
#include "dc.h"
#include <iostream>
#include <cmath>

using namespace std;

dc::dc(const circuit * const c) : analysis(c) {

    // count number of each element
    size_t n = c->nodes.size();
    size_t m = 0;
    for (const auto & e : c->linelems) {
        switch (e->ElemType) {
            case circuit::linelem::C:
            case circuit::linelem::V:
                if (e->Node1==e->Node2) break;
                if (std::isnan(voltage(e))) break;
                V_source_i[e] = m;
                m++;
            default: break;
        }
    }
    for (const auto & e : itrelems) {
        switch (e->ElemType) {
            case circuit::linelem::C:
            case circuit::linelem::V:
                if (e->Node1==e->Node2) break;
                if (std::isnan(voltage(e))) break;
                V_source_i[e] = m;
                m++;
            default: break;
        }
    }

    // A and z matricies
    Eigen::SparseMatrix<double> A(n+m, n+m);
    Eigen::SparseMatrix<double> z(n+m, 1);
    for (const auto & e : c->linelems) {
        stamp_A(A, e);
        stamp_z(z, e);
    }
    for (const auto & e : itrelems) {
        stamp_A(A, e);
        stamp_z(z, e);
    }

    Eigen::SparseLU< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > A_solver;
    Eigen::SparseMatrix<double> x(n+m, 1);
    struct mosfet_nodes{double vd, vs, vg;};
    unordered_map<circuit::mosfet*,mosfet_nodes> k;
    for (auto m : c->mosfets) k[m] = {0, 0, 0};
    for (size_t nr_i = 0; nr_i < 100; nr_i++) {
        Eigen::SparseMatrix<double> NR_A(n+m, n+m);
        Eigen::SparseMatrix<double> NR_z(n+m, 1);
        for (auto m : c->mosfets) {
            double V_GS = k[m].vg-k[m].vs;
            double V_DS = k[m].vd-k[m].vs;
            stamp_A_NR(NR_A, m, V_GS, V_DS);
            stamp_z_NR(NR_z, m, V_GS, V_DS);
        }
        Eigen::SparseMatrix<double> A_combined = A + NR_A;
        A_solver.analyzePattern(A_combined);
        A_solver.factorize(A_combined);
        while (A_solver.info()) {
            cerr << "Could not factorize A in DC." << endl
                << Eigen::MatrixXd(A) << endl
                << Eigen::MatrixXd(NR_A) << endl
                << Eigen::MatrixXd(A_combined) << endl
                << " Exiting..." << endl;
            exit(1);
        }

        Eigen::SparseMatrix<double> z_combined = z + NR_z;
        x = A_solver.solve(z_combined);

        bool changed = false;
        for (auto m : c->mosfets) {
            double vs = k[m].vs;
            double vd = k[m].vd;
            double vg = k[m].vg;
            k[m].vs = ((m->NodeS==circuit::gnd) ? 0 : x.coeffRef(m->NodeS->i, 0));
            k[m].vd = ((m->NodeD==circuit::gnd) ? 0 : x.coeffRef(m->NodeD->i, 0));
            k[m].vg = ((m->NodeG==circuit::gnd) ? 0 : x.coeffRef(m->NodeG->i, 0));
            changed = changed || (abs(vs-k[m].vs)>=precision) || (abs(vd-k[m].vd)>=precision) || (abs(vg-k[m].vg)>=precision);
        }
        if (!changed) {
            break;
        }
    }

    // record voltages
    for (auto node_i : c->nodes) {
        node_voltage[node_i.second] = x.coeffRef(node_i.second->i,0 );
    }

    // record currents
    for (auto e : c->linelems) {
        switch (e->ElemType) {
            case circuit::linelem::C: {
                auto e_c = (circuit::capacitor*)e;
                if (V_source_i.find(e)==V_source_i.end()) capacitor_current[e_c] = 0;
                capacitor_current[e_c] = x.coeffRef(n+V_source_i.at(e_c), 0);
                } break;
            case circuit::linelem::V: {
                auto e_v = (circuit::V_source*)e;
                if (V_source_i.find(e)==V_source_i.end()) V_source_current[e_v] = 0;
                V_source_current[e_v] = x.coeffRef(n+V_source_i.at(e_v), 0);
                } break;
            default: break;
        }
    }
    for (auto e : itrelems) {
        switch (e->ElemType) {
            case circuit::linelem::C: {
                auto e_c = (circuit::capacitor*)e;
                if (V_source_i.find(e)==V_source_i.end()) capacitor_current[e_c] = 0;
                else capacitor_current[e_c] = x.coeffRef(n+V_source_i.at(e_c), 0);
                } break;
            case circuit::linelem::V: {
                auto e_v = (circuit::V_source*)e;
                if (V_source_i.find(e)==V_source_i.end()) V_source_current[e_v] = 0;
                else V_source_current[e_v] = x.coeffRef(n+V_source_i.at(e_v), 0);
                } break;
            default: break;
        }
    }

}


double dc::mosfet_solver::conductance(const circuit::mosfet * const m, const double & V_GS, const double & V_DS) {
    if (m->ElemType==circuit::mosfet::nmos) {
        // cutoff
        if (V_GS <= m->V_T)
            return 0;

        // extra case for correct convergence
        if (V_DS < 0)
            return 1;

        double V_GST = V_GS-m->V_T;
        double k = m->MU * m->C_OX * m->W / m->L;
        double gm = 1 + m->LAMBDA * V_DS;

        // linear
        if (V_DS <= V_GST)
            return k*m->LAMBDA*(V_GST*V_DS-pow(V_DS,2)/2) + k*(-V_DS+V_GST)*gm;

        // saturation
        return 0.5 * k * pow(V_GST,2) * m->LAMBDA;
    } else if (m->ElemType==circuit::mosfet::pmos) {
        // cutoff
        if (V_GS >= m->V_T)
            return 0;

        // extra case for correct convergence
        if (V_DS > 0)
            return 1;

        double V_GST = V_GS-m->V_T;
        double k = m->MU * m->C_OX * m->W / m->L;
        double gm = 1 - m->LAMBDA * V_DS;

        // linear
        if (V_DS >= V_GST)
            return k*m->LAMBDA*(V_GST*V_DS-pow(V_DS,2)/2) - k*(V_GST-V_DS)*gm;

        // saturation
        return 0.5 * k * pow(V_GST,2) * m->LAMBDA;
    } else {
        cerr << "Unknown MOSFET type: " << m->ElemType << endl;
        throw circuit::mosfet::UnsupportedMOSFETType();
    }
}
double dc::mosfet_solver::I_DS(const circuit::mosfet * const m, const double & V_GS, const double & V_DS) {
    if (m->ElemType==circuit::mosfet::nmos) {
        // cutoff
        if (V_GS <= m->V_T)
            return 0;

        // extra case for correct convergence
        if (V_DS < 0)
            return V_DS;

        double V_GST = V_GS-m->V_T;
        double k = m->MU * m->C_OX * m->W / m->L;
        double gm = 1 + (m->LAMBDA*V_DS);

        // linear
        if (V_DS <= V_GST) {
            return k * (V_DS*V_GST - 0.5*pow(V_DS,2)) * gm;
        }

        // saturation
        return 0.5 * k * pow(V_GST,2) * gm;
    } else if (m->ElemType==circuit::mosfet::pmos) {
        // cutoff
        if (V_GS >= m->V_T)
            return 0;

        // extra case for correct convergence
        if (V_DS > 0)
            return V_DS;

        double V_GST = V_GS-m->V_T;
        double k = m->MU * m->C_OX * m->W / m->L;
        double gm = 1 - m->LAMBDA * V_DS;

        // linear
        if (V_DS >= V_GST)
            return -k * (V_DS*V_GST - 0.5*pow(V_DS,2)) * gm;

        // saturation
        return -0.5 * k * pow(V_GST,2) * gm;
    } else {
        cerr << "Unknown MOSFET type: " << m->ElemType << endl;
        throw circuit::mosfet::UnsupportedMOSFETType();
    }
}

double dc::mosfet_solver::NR_G_eq(const circuit::mosfet * const m, const double & V_GS,const double & V_DS) {
    return mosfet_solver::conductance(m, V_GS, V_DS);
}

double dc::mosfet_solver::NR_I_eq(const circuit::mosfet * const m, const double & V_GS,const double & V_DS) {
    return mosfet_solver::I_DS(m, V_GS, V_DS) - (mosfet_solver::conductance(m, V_GS, V_DS)*V_DS);
}


void dc::stamp_A(Eigen::SparseMatrix<double> & A, const circuit::linelem * const e) const {
    switch (e->ElemType) {
        case circuit::linelem::R: {
            auto e_r = (circuit::resistor*)e;
            if (e_r->Node1!=circuit::gnd)
                A.coeffRef(e_r->Node1->i, e_r->Node1->i) += 1.0 / e_r->resistance;
            if (e_r->Node2!=circuit::gnd)
                A.coeffRef(e_r->Node2->i, e_r->Node2->i) += 1.0 / e_r->resistance;
            if (e_r->Node1!=circuit::gnd && e_r->Node2!=circuit::gnd) {
                A.coeffRef(e_r->Node1->i, e_r->Node2->i) -= 1.0 / e_r->resistance;
                A.coeffRef(e_r->Node2->i, e_r->Node1->i) -= 1.0 / e_r->resistance;
            }
            } break;
        case circuit::linelem::V: {
            const auto e_v = (circuit::V_source*)e;
            if (V_source_i.find(e)==V_source_i.end()) break;
            if (e_v->Node1!=circuit::gnd) {
                A.coeffRef(e_v->Node1->i, c->nodes.size()+V_source_i.at(e)) = -1.0;
                A.coeffRef(c->nodes.size()+V_source_i.at(e), e_v->Node1->i) = -1.0;
            }
            if (e_v->Node2!=circuit::gnd) {
                A.coeffRef(e_v->Node2->i, c->nodes.size()+V_source_i.at(e)) = 1.0;
                A.coeffRef(c->nodes.size()+V_source_i.at(e), e_v->Node2->i) = 1.0;
            }
            } break;
        case circuit::linelem::C: {
            const auto e_c = (circuit::capacitor*)e;
            if (V_source_i.find(e)==V_source_i.end()) break;
            if (e_c->Node1!=circuit::gnd) {
                A.coeffRef(e_c->Node1->i, c->nodes.size()+V_source_i.at(e)) = -1.0;
                A.coeffRef(c->nodes.size()+V_source_i.at(e), e_c->Node1->i) = -1.0;
            }
            if (e_c->Node2!=circuit::gnd) {
                A.coeffRef(e_c->Node2->i, c->nodes.size()+V_source_i.at(e)) = 1.0;
                A.coeffRef(c->nodes.size()+V_source_i.at(e), e_c->Node2->i) = 1.0;
            }
            } break;
        default: break;
    }
}

void dc::stamp_z(Eigen::SparseMatrix<double> & z, const circuit::linelem * const e) const {
    switch (e->ElemType) {
        case circuit::linelem::C: {
            auto e_c = (circuit::capacitor*)e;
            if (V_source_i.find(e)==V_source_i.end()) break;
            z.coeffRef(c->nodes.size()+V_source_i.at(e), 0) = -1.0 * voltage(e_c);
            } break;
        case circuit::linelem::L: {
            auto e_l = (circuit::inductor*)e;
            if (e_l->Node1!=circuit::gnd)
                z.coeffRef(e_l->Node1->i, 0) += current(e_l);
            if (e_l->Node2!=circuit::gnd)
                z.coeffRef(e_l->Node2->i, 0) -= current(e_l);
            } break;
        case circuit::linelem::V: {
            auto e_v = (circuit::V_source*)e;
            if (V_source_i.find(e)==V_source_i.end()) break;
            z.coeffRef(c->nodes.size()+V_source_i.at(e), 0) = -1.0 * voltage(e_v);
            } break;
        case circuit::linelem::I: {
            auto e_i = (circuit::I_source*)e;
            if (e_i->Node1!=circuit::gnd)
                z.coeffRef(e_i->Node1->i, 0) += current(e_i);
            if (e_i->Node2!=circuit::gnd)
                z.coeffRef(e_i->Node2->i, 0) -= current(e_i);
            } break;
        default: break;
    }
}
void dc::stamp_A_NR(Eigen::SparseMatrix<double> & A_NR, const circuit::mosfet * const m, const double & V_GS, const double & V_DS) const {
    if (m->NodeS!=circuit::gnd) {
        A_NR.coeffRef(m->NodeS->i, m->NodeS->i) += mosfet_solver::NR_G_eq(m, V_GS, V_DS);
    }
    if (m->NodeD!=circuit::gnd) {
        A_NR.coeffRef(m->NodeD->i, m->NodeD->i) += mosfet_solver::NR_G_eq(m, V_GS, V_DS);
    }
    if (m->NodeS!=circuit::gnd && m->NodeD!=circuit::gnd) {
        A_NR.coeffRef(m->NodeD->i, m->NodeS->i) -= mosfet_solver::NR_G_eq(m, V_GS, V_DS);
        A_NR.coeffRef(m->NodeS->i, m->NodeD->i) -= mosfet_solver::NR_G_eq(m, V_GS, V_DS);
    }
}
void dc::stamp_z_NR(Eigen::SparseMatrix<double> & z_NR, const circuit::mosfet * const m, const double & V_GS, const double & V_DS) const {
    if (m->NodeS!=circuit::gnd) {
        z_NR.coeffRef(m->NodeS->i, 0) += mosfet_solver::NR_I_eq(m, V_GS, V_DS);
    }
    if (m->NodeD!=circuit::gnd) {
        z_NR.coeffRef(m->NodeD->i, 0) -= mosfet_solver::NR_I_eq(m, V_GS, V_DS);
    }
}


void dc::plotnv(matlab * const m, const int & node_name) const {
    std::cerr << "Warning: Cannot plot DC voltage" << endl;
}

void dc::printnv(const int & node_name) const {
    circuit::node * n = 0;
    for (const auto i : c->nodes)
        if (i.second->name == node_name)
            n = i.second;
    if (n==0) {
        std::cerr << "Warning: Could not find node " << node_name << endl;
        return;
    }
    cout << "Node " << node_name << " voltage: " << node_voltage.at(n) << endl;
}
