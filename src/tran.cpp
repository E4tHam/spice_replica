
#include "tran.h"
#include <iostream>

using namespace std;

tran::tran(const circuit * const c, const double & time_step, const double & stop_time)
    : analysis(c), time_step(time_step) {

    // set n
    size_t n = c->nodes.size();

    // set voltage source indicies for MNA, and set m
    size_t m = 0;
    for (const auto & e : c->linelems) {
        if (e->ElemType==circuit::linelem::V) {
            V_source_i[(circuit::V_source*)e] = m;
            m++;
        }
    }
    for (const auto & e : itrelems) {
        if (e->ElemType==circuit::linelem::V) {
            V_source_i[(circuit::V_source*)e] = m;
            m++;
        }
    }

    // dc analysis for initial state
    dc initial_dc(c);
    for (const auto & v : initial_dc.node_voltage) {
        node_voltage[v.first].push_back(v.second);
    }
    for (const auto & i : initial_dc.capacitor_current) {
        storage_device_current[(circuit::storage_device*)i.first].push_back(i.second);
    }
    for (const auto & e : c->linelems) if (e->ElemType==circuit::linelem::L) {
        storage_device_current[(circuit::storage_device*)e].push_back(initial_dc.current(e));
    }
    for (const auto & e : itrelems) if (e->ElemType==circuit::linelem::L) {
        storage_device_current[(circuit::storage_device*)e].push_back(initial_dc.current(e));
    }
    for (const auto & i : initial_dc.V_source_current)
        V_source_current[i.first].push_back(i.second);

    // Create A matrix
    Eigen::SparseMatrix<double> A(n+m, n+m);
    for (const auto & e : c->linelems) {
        stamp_A(A, e);
    }
    for (const auto & e : itrelems) {
        stamp_A(A, e);
    }

    // Run each time step
    this->stop_time = 0;
    while (this->stop_time < stop_time) {
        // Create z matrix
        Eigen::SparseMatrix<double> z(n+m, 1);
        for (const auto & e : c->linelems) {
            stamp_z(z, e);
        }
        for (const auto & e : itrelems) {
            stamp_z(z, e);
        }

        // solve x matrix
        Eigen::SparseLU< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > A_solver;
        Eigen::SparseMatrix<double> x(n+m, 1);
        if (c->mosfets.size()==0) {
            A_solver.analyzePattern(A);
            A_solver.factorize(A);
            x = A_solver.solve(z);
        } else {
            // todo: NR
            struct mosfet_nodes{double vd, vs, vg;};
            unordered_map<circuit::mosfet*, mosfet_nodes> k;
            for (const auto & m : c->mosfets) k[m] = {voltage(m->NodeD), voltage(m->NodeS), voltage(m->NodeG)};
            for (size_t nr_i = 0; nr_i < 100; nr_i++) {
                Eigen::SparseMatrix<double> NR_A(n+m, n+m);
                Eigen::SparseMatrix<double> NR_z(n+m, 1);
                for (const auto & m : c->mosfets) {
                    double V_GS = k[m].vg-k[m].vs;
                    double V_DS = k[m].vd-k[m].vs;
                    stamp_A_NR(NR_A, m, V_GS, V_DS);
                    stamp_z_NR(NR_z, m, V_GS, V_DS);
                }
                Eigen::SparseMatrix<double> A_combined = A + NR_A;
                A_solver.analyzePattern(A_combined);
                A_solver.factorize(A_combined);
                while (A_solver.info()) {
                    cerr << "Could not factorize A in TRAN." << endl << Eigen::MatrixXd(A_combined) << endl << " Exiting..." << endl;
                    exit(1);
                }

                Eigen::SparseMatrix<double> z_combined = z + NR_z;
                x = A_solver.solve(z_combined);

                bool changed = false;
                for (const auto & m : c->mosfets) {
                    double vs = k[m].vs;
                    double vd = k[m].vd;
                    double vg = k[m].vg;
                    k[m].vs = ((m->NodeS==circuit::gnd) ? 0 : x.coeffRef(m->NodeS->i, 0));
                    k[m].vd = ((m->NodeD==circuit::gnd) ? 0 : x.coeffRef(m->NodeD->i, 0));
                    k[m].vg = ((m->NodeG==circuit::gnd) ? 0 : x.coeffRef(m->NodeG->i, 0));
                    changed = changed || (abs(vs-k[m].vs)>precision) || (abs(vd-k[m].vd)>precision) || (abs(vg-k[m].vg)>precision);
                }
                if (!changed) {
                    break;
                }
            }
        }

        // record voltages
        for (const auto & n : c->nodes) {
            node_voltage[n.second].push_back( x.coeffRef(n.second->i,0 ) );
        }

        // record currents (could be clearer)
        for (const auto & e : c->linelems) {
            extract_current(x, e);
        }
        for (const auto & e : itrelems) {
            extract_current(x, e);
        }

        this->stop_time += time_step;
    }

}

void tran::plotnv(matlab * const m, const int & node_name) const {
    circuit::node * n = 0;
    for (const auto i : c->nodes)
        if (i.second->name == node_name)
            n = i.second;
    if (n==0) {
        std::cerr << "Warning: Could not find node " << node_name << endl;
        return;
    }
    m->show_plot(node_voltage.at(n), ("Node "+std::to_string(node_name)+" Voltage"), "time (s)", "Voltage", time_step, stop_time);
}

void tran::printnv(const int & node_name) const {
    circuit::node * n = 0;
    for (const auto i : c->nodes)
        if (i.second->name == node_name)
            n = i.second;
    if (n==0) {
        std::cerr << "Warning: Could not find node " << node_name << endl;
        return;
    }
    cout << "Node " << node_name << " voltages: ";
    for (const auto v : node_voltage.at(n)) {
        cout << v << " ";
    }
    cout << endl;
}


void tran::stamp_A(Eigen::SparseMatrix<double> & A, const circuit::linelem * const e) const {
    switch (e->ElemType) {
        case circuit::linelem::R: {
            const auto e_r = (circuit::resistor*)e;
            if (e_r->Node1!=circuit::gnd)
                A.coeffRef(e_r->Node1->i, e_r->Node1->i) += 1.0 / e_r->resistance;
            if (e_r->Node2!=circuit::gnd)
                A.coeffRef(e_r->Node2->i, e_r->Node2->i) += 1.0 / e_r->resistance;
            if (e_r->Node1!=circuit::gnd && e_r->Node2!=circuit::gnd) {
                A.coeffRef(e_r->Node1->i, e_r->Node2->i) -= 1.0 / e_r->resistance;
                A.coeffRef(e_r->Node2->i, e_r->Node1->i) -= 1.0 / e_r->resistance;
            }
        } break;
        case circuit::linelem::C: {
            const auto e_c = (circuit::capacitor*)e;
            if (e_c->Node1!=circuit::gnd) {
                A.coeffRef(e_c->Node1->i, e_c->Node1->i) += storage_device_solver::conductance(e_c, time_step);
            }
            if (e_c->Node2!=circuit::gnd) {
                A.coeffRef(e_c->Node2->i, e_c->Node2->i) += storage_device_solver::conductance(e_c, time_step);
            }
            if (e_c->Node1!=circuit::gnd && e_c->Node2!=circuit::gnd) {
                A.coeffRef(e_c->Node1->i, e_c->Node2->i) -= storage_device_solver::conductance(e_c, time_step);
                A.coeffRef(e_c->Node2->i, e_c->Node1->i) -= storage_device_solver::conductance(e_c, time_step);
            }
            } break;
        case circuit::linelem::L: {
            const auto e_l = (circuit::inductor*)e;
            if (e_l->Node1!=circuit::gnd) {
                A.coeffRef(e_l->Node1->i, e_l->Node1->i) += storage_device_solver::conductance(e_l, time_step);
            }
            if (e_l->Node2!=circuit::gnd) {
                A.coeffRef(e_l->Node2->i, e_l->Node2->i) += storage_device_solver::conductance(e_l, time_step);
            }
            if (e_l->Node1!=circuit::gnd && e_l->Node2!=circuit::gnd) {
                A.coeffRef(e_l->Node1->i, e_l->Node2->i) -= storage_device_solver::conductance(e_l, time_step);
                A.coeffRef(e_l->Node2->i, e_l->Node1->i) -= storage_device_solver::conductance(e_l, time_step);
            }
            } break;
        case circuit::linelem::V: {
            const auto e_v = (circuit::V_source*)e;
            if (e_v->Node1!=circuit::gnd) {
                A.coeffRef(e_v->Node1->i, c->nodes.size()+V_source_i.at(e_v)) = -1.0;
                A.coeffRef(c->nodes.size()+V_source_i.at(e_v), e_v->Node1->i) = -1.0;
            }
            if (e_v->Node2!=circuit::gnd) {
                A.coeffRef(e_v->Node2->i, c->nodes.size()+V_source_i.at(e_v)) = 1.0;
                A.coeffRef(c->nodes.size()+V_source_i.at(e_v), e_v->Node2->i) = 1.0;
            }
            } break;
        default: break;
    }
}


// run after recording all node voltages
void tran::extract_current(Eigen::SparseMatrix<double> & x, const circuit::linelem * const e) {
    switch (e->ElemType) {
        case circuit::linelem::C: {
            auto e_c = (circuit::capacitor*)e;
            double I = (storage_device_solver::conductance(e_c, time_step)*voltage(e_c, -2)) - (storage_device_solver::conductance(e_c, time_step)*voltage(e_c, -1));
            storage_device_current[e_c].push_back( I );
            } break;
        case circuit::linelem::L: {
            auto e_l = (circuit::inductor*)e;
            double I = current(e_l) - (storage_device_solver::conductance(e_l, time_step) * voltage(e_l, -1));
            storage_device_current[e_l].push_back( I );
            } break;
        case circuit::linelem::V: {
            auto e_v = (circuit::V_source*)e;
            double I = x.coeffRef(c->nodes.size()+V_source_i.at(e_v), 0);
            V_source_current[e_v].push_back( I );
            } break;
        default: break;
    }
}



void tran::stamp_z(Eigen::SparseMatrix<double> & z, const circuit::linelem * const e) const {
    switch (e->ElemType) {
        case circuit::linelem::C: {
            auto e_c = (circuit::capacitor*)e;
            double I_eq = storage_device_solver::conductance(e_c, time_step) * voltage(e_c);
            if (e_c->Node1!=circuit::gnd)
                z.coeffRef(e_c->Node1->i, 0) += I_eq;
            if (e_c->Node2!=circuit::gnd)
                z.coeffRef(e_c->Node2->i, 0) -= I_eq;
            } break;
        case circuit::linelem::L: {
            auto e_l = (circuit::inductor*)e;
            double I_eq = current(e_l);
            if (e_l->Node1!=circuit::gnd)
                z.coeffRef(e_l->Node1->i, 0) += I_eq;
            if (e_l->Node2!=circuit::gnd)
                z.coeffRef(e_l->Node2->i, 0) -= I_eq;
            } break;
        case circuit::linelem::V: {
            auto e_v = (circuit::V_source*)e;
            z.coeffRef(c->nodes.size()+V_source_i.at(e_v), 0) = -1.0 * voltage(e_v);
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




void tran::stamp_A_NR(Eigen::SparseMatrix<double> & A_NR, const circuit::mosfet * const m, const double & V_GS, const double & V_DS) const {
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
void tran::stamp_z_NR(Eigen::SparseMatrix<double> & z_NR, const circuit::mosfet * const m, const double & V_GS, const double & V_DS) const {
    if (m->NodeS!=circuit::gnd) {
        z_NR.coeffRef(m->NodeS->i, 0) += mosfet_solver::NR_I_eq(m, V_GS, V_DS);
    }
    if (m->NodeD!=circuit::gnd) {
        z_NR.coeffRef(m->NodeD->i, 0) -= mosfet_solver::NR_I_eq(m, V_GS, V_DS);
    }
}





// MOSFET Equations

double tran::mosfet_solver::conductance(const circuit::mosfet * const m, const double & V_GS, const double & V_DS) {
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
double tran::mosfet_solver::I_DS(const circuit::mosfet * const m, const double & V_GS, const double & V_DS) {
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

double tran::mosfet_solver::NR_G_eq(const circuit::mosfet * const m, const double & V_GS,const double & V_DS) {
    return mosfet_solver::conductance(m, V_GS, V_DS);
}

double tran::mosfet_solver::NR_I_eq(const circuit::mosfet * const m, const double & V_GS,const double & V_DS) {
    return mosfet_solver::I_DS(m, V_GS, V_DS) - (mosfet_solver::conductance(m, V_GS, V_DS)*V_DS);
}
