
#include "dc.h"
#include <iostream>

using namespace std;

dc::dc(const circuit * const c) : analysis(c) {

    // count number of each element
    size_t n, m;
    size_t num_capacitors = 0;
    size_t num_voltage_sources = 0;
    for (auto e : c->linelems) {
        switch (e->ElemType) {
            case circuit::linelem::C: num_capacitors++; break;
            case circuit::linelem::V: num_voltage_sources++; break;
            default: break;
        }
    }

    // temporary matricies
    Eigen::MatrixXd A_dense, B, C, D, G, x, NR_G;
    Eigen::Matrix <double, Eigen::Dynamic, 1> I, E, z, NR_I;
    size_t E_i;

    // set G array with resistors only
    n = c->nodes.size();
    G = Eigen::MatrixXd::Zero(n, n);
    for (auto e : c->linelems) {
        if (e->ElemType == circuit::linelem::R) {
            auto e_r = (circuit::resistor*)e;
            if (e_r->Node1!=circuit::gnd)
                G(e_r->Node1->i, e_r->Node1->i) += 1.0 / e_r->resistance;
            if (e_r->Node2!=circuit::gnd)
                G(e_r->Node2->i, e_r->Node2->i) += 1.0 / e_r->resistance;
            if (e_r->Node1!=circuit::gnd && e_r->Node2!=circuit::gnd) {
                G(e_r->Node1->i, e_r->Node2->i) -= 1.0 / e_r->resistance;
                G(e_r->Node2->i, e_r->Node1->i) -= 1.0 / e_r->resistance;
            }
        }
    }

    // Find remaining matricies modeling storage elements as power sources
    m = num_capacitors + num_voltage_sources;
    B = Eigen::MatrixXd::Zero(n, m);
    I = Eigen::MatrixXd::Zero(n, 1);
    E = Eigen::MatrixXd::Zero(m, 1);
    E_i = 0;
    for (auto e : c->linelems) {
        switch (e->ElemType) {
            case circuit::linelem::C: {
                auto e_c = (circuit::capacitor*)e;
                if (e_c->Node1!=circuit::gnd)
                    B(e_c->Node1->i, E_i) = -1.0;
                if (e_c->Node2!=circuit::gnd)
                    B(e_c->Node2->i, E_i) = 1.0;
                E(E_i, 0) = -1.0 * e_c->initial_voltage;
                E_i++;
                } break;
            case circuit::linelem::L: {
                auto e_l = (circuit::inductor*)e;
                if (e_l->Node1!=circuit::gnd)
                    I(e_l->Node1->i, 0) += e_l->initial_current;
                if (e_l->Node2!=circuit::gnd)
                    I(e_l->Node2->i, 0) -= e_l->initial_current;
                } break;
            case circuit::linelem::V: {
                auto e_v = (circuit::V_source*)e;
                if (e_v->Node1!=circuit::gnd)
                    B(e_v->Node1->i, E_i) = -1.0;
                if (e_v->Node2!=circuit::gnd)
                    B(e_v->Node2->i, E_i) = 1.0;
                E(E_i, 0) = -1.0 * voltage(e_v);
                E_i++;
                } break;
            case circuit::linelem::I: {
                auto e_i = (circuit::I_source*)e;
                if (e_i->Node1!=circuit::gnd)
                    I(e_i->Node1->i, 0) += current(e_i);
                if (e_i->Node2!=circuit::gnd)
                    I(e_i->Node2->i, 0) -= current(e_i);
                } break;
            default: break;
        }
    }
    C = B.transpose();
    D = Eigen::MatrixXd::Zero(m, m);
    A_dense = Eigen::MatrixXd(n+m, n+m);
    A_dense << G, B, C, D;
    z = Eigen::MatrixXd(n+m, 1);

    struct mosfet_nodes{double vd, vs, vg;};
    unordered_map<circuit::mosfet*,mosfet_nodes> k;
    for (auto m : c->mosfets) k[m] = {0, 0, 0};
    for (size_t nr_i = 0; nr_i < 100; nr_i++) {
        NR_G = Eigen::MatrixXd::Zero(n, n);
        NR_I = Eigen::MatrixXd::Zero(n, 1);
        for (auto m : c->mosfets) {
            double V_GS = k[m].vg-k[m].vs;
            double V_DS = k[m].vd-k[m].vs;
            if (m->NodeS!=circuit::gnd) {
                NR_I(m->NodeS->i, 0) += mosfet_solver::NR_I_eq(m, V_GS, V_DS);
                NR_G(m->NodeS->i, m->NodeS->i) += mosfet_solver::NR_G_eq(m, V_GS, V_DS);
            }
            if (m->NodeD!=circuit::gnd) {
                NR_I(m->NodeD->i, 0) -= mosfet_solver::NR_I_eq(m, V_GS, V_DS);
                NR_G(m->NodeD->i, m->NodeD->i) += mosfet_solver::NR_G_eq(m, V_GS, V_DS);
            }
            if (m->NodeS!=circuit::gnd && m->NodeD!=circuit::gnd) {
                NR_G(m->NodeD->i, m->NodeS->i) -= mosfet_solver::NR_G_eq(m, V_GS, V_DS);
                NR_G(m->NodeS->i, m->NodeD->i) -= mosfet_solver::NR_G_eq(m, V_GS, V_DS);
            }
        }
        A_dense << (G+NR_G), B, C, D;
        z << (I+NR_I), E;
        x = A_dense.completeOrthogonalDecomposition().solve(z);

        bool changed = false;
        for (auto m : c->mosfets) {
            double vs = k[m].vs;
            double vd = k[m].vd;
            double vg = k[m].vg;
            k[m].vs = ((m->NodeS==circuit::gnd) ? 0 : x(m->NodeS->i, 0));
            k[m].vd = ((m->NodeD==circuit::gnd) ? 0 : x(m->NodeD->i, 0));
            k[m].vg = ((m->NodeG==circuit::gnd) ? 0 : x(m->NodeG->i, 0));
            changed = changed || (abs(vs-k[m].vs)>=precision) || (abs(vd-k[m].vd)>=precision) || (abs(vg-k[m].vg)>=precision);
        }
        if (!changed) {
            break;
        }
    }

    // record voltages
    for (auto node_i : c->nodes) {
        node_voltage[node_i.second] = x(node_i.second->i,0 );
    }

    // record currents
    E_i = 0;
    for (auto e : c->linelems) {
        switch (e->ElemType) {
            case circuit::linelem::C: {
                auto e_c = (circuit::capacitor*)e;
                capacitor_current[e_c] = x(n+E_i, 0);
                E_i++;
                } break;
            case circuit::linelem::V: {
                auto e_v = (circuit::V_source*)e;
                V_source_current[e_v] = x(n+E_i, 0);
                E_i++;
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
        exit(1);
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
        exit(1);
    }
}

double dc::mosfet_solver::NR_G_eq(const circuit::mosfet * const m, const double & V_GS,const double & V_DS) {
    return mosfet_solver::conductance(m, V_GS, V_DS);
}

double dc::mosfet_solver::NR_I_eq(const circuit::mosfet * const m, const double & V_GS,const double & V_DS) {
    return mosfet_solver::I_DS(m, V_GS, V_DS) - (mosfet_solver::conductance(m, V_GS, V_DS)*V_DS);
}
