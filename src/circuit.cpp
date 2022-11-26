
#include "circuit.h"
#include "circuit_interface.h"
#include <iostream>

using namespace std;



// static variables
circuit::node GND_INSTANCE(-1, 0, -1);
circuit::node * const circuit::gnd = &GND_INSTANCE;





// circuit methods
circuit::circuit() { }

circuit::circuit(const std::string & js) {

    // initialize this->nodes and this->linelems
    circuit_interface::circuit_from_json(this, js);

    tran_precision = 0.000000001;

}







circuit::~circuit() {
    for (auto e : linelems)
        if (e) delete e;
    for (auto n : nodes)
        if (n.second) delete n.second;
}




void circuit::dc() {

    step_num = 0;

    // count number of each element
    size_t n, m;
    size_t num_capacitors = 0;
    size_t num_voltage_sources = 0;
    for (auto e : linelems) {
        switch (e->ElemType) {
            case linelem::C: num_capacitors++; break;
            case linelem::V: num_voltage_sources++; break;
            default: break;
        }
    }

    // temporary matricies
    Eigen::MatrixXd A_dense, B, C, D, G, x, NR_G;
    Eigen::Matrix <double, Eigen::Dynamic, 1> I, E, z, NR_I;
    size_t E_i;

    // set G array with resistors only
    n = nodes.size();
    G = Eigen::MatrixXd::Zero(n, n);
    for (auto e : linelems) {
        if (e->ElemType == linelem::R) {
            auto e_r = (resistor*)e;
            if (e_r->Node1!=gnd)
                G(e_r->Node1->i, e_r->Node1->i) += 1.0 / e_r->resistance;
            if (e_r->Node2!=gnd)
                G(e_r->Node2->i, e_r->Node2->i) += 1.0 / e_r->resistance;
            if (e_r->Node1!=gnd && e_r->Node2!=gnd) {
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
    for (auto e : linelems) {
        switch (e->ElemType) {
            case linelem::C: {
                auto e_c = (capacitor*)e;
                if (e_c->Node1!=gnd)
                    B(e_c->Node1->i, E_i) = -1.0;
                if (e_c->Node2!=gnd)
                    B(e_c->Node2->i, E_i) = 1.0;
                E(E_i, 0) = -1.0 * e_c->initial_voltage;
                E_i++;
                } break;
            case linelem::L: {
                auto e_l = (inductor*)e;
                if (e_l->Node1!=gnd)
                    I(e_l->Node1->i, 0) += e_l->initial_current;
                if (e_l->Node2!=gnd)
                    I(e_l->Node2->i, 0) -= e_l->initial_current;
                } break;
            case linelem::V: {
                auto e_v = (V_source*)e;
                if (e_v->Node1!=gnd)
                    B(e_v->Node1->i, E_i) = -1.0;
                if (e_v->Node2!=gnd)
                    B(e_v->Node2->i, E_i) = 1.0;
                E(E_i, 0) = -1.0 * e_v->voltage();
                E_i++;
                } break;
            case linelem::I: {
                auto e_i = (I_source*)e;
                if (e_i->Node1!=gnd)
                    I(e_i->Node1->i, 0) += e_i->current();
                if (e_i->Node2!=gnd)
                    I(e_i->Node2->i, 0) -= e_i->current();
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
    unordered_map<mosfet*,mosfet_nodes> k;
    for (auto m : mosfets) k[m] = {0, 0, 0};
    for (size_t nr_i = 0; nr_i < 100; nr_i++) {
        NR_G = Eigen::MatrixXd::Zero(n, n);
        NR_I = Eigen::MatrixXd::Zero(n, 1);
        for (auto m : mosfets) {
            double V_GS = k[m].vg-k[m].vs;
            double V_DS = k[m].vd-k[m].vs;
            if (m->NodeS!=gnd) {
                NR_I(m->NodeS->i, 0) += m->NR_I_eq(V_GS, V_DS);
                NR_G(m->NodeS->i, m->NodeS->i) += m->NR_G_eq(V_GS, V_DS);
            }
            if (m->NodeD!=gnd) {
                NR_I(m->NodeD->i, 0) -= m->NR_I_eq(V_GS, V_DS);
                NR_G(m->NodeD->i, m->NodeD->i) += m->NR_G_eq(V_GS, V_DS);
            }
            if (m->NodeS!=gnd && m->NodeD!=gnd) {
                NR_G(m->NodeD->i, m->NodeS->i) -= m->NR_G_eq(V_GS, V_DS);
                NR_G(m->NodeS->i, m->NodeD->i) -= m->NR_G_eq(V_GS, V_DS);
            }
        }
        A_dense << (G+NR_G), B, C, D;
        z << (I+NR_I), E;
        x = A_dense.completeOrthogonalDecomposition().solve(z);

        bool changed = false;
        for (auto m : mosfets) {
            double vs = k[m].vs;
            double vd = k[m].vd;
            double vg = k[m].vg;
            k[m].vs = ((m->NodeS==gnd) ? 0 : x(m->NodeS->i, 0));
            k[m].vd = ((m->NodeD==gnd) ? 0 : x(m->NodeD->i, 0));
            k[m].vg = ((m->NodeG==gnd) ? 0 : x(m->NodeG->i, 0));
            changed = changed || (abs(vs-k[m].vs)>=tran_precision) || (abs(vd-k[m].vd)>=tran_precision) || (abs(vg-k[m].vg)>=tran_precision);
        }
        if (!changed) {
            // cout << "exited tran at nr_i=" << nr_i << endl;
            break;
        }
    }

    // cout << A_dense << endl << x << endl << z << endl;

    // record voltages
    for (auto node_i : nodes) {
        node_i.second->voltages.push_back( x(node_i.second->i,0 ) );
    }

    // record currents
    E_i = 0;
    for (auto e : linelems) {
        switch (e->ElemType) {
            case linelem::C: {
                auto e_c = (capacitor*)e;
                e_c->currents.push_back( x(n+E_i, 0) );
                E_i++;
                } break;
            case linelem::L: {
                auto e_l = (inductor*)e;
                e_l->currents.push_back( e_l->initial_current );
                } break;
            case linelem::V: {
                auto e_v = (V_source*)e;
                e_v->currents.push_back( x(n+E_i, 0) );
                E_i++;
                } break;
            default: break;
        }
    }
    // for (auto e : itrelems) {
    //     switch (e->ElemType) {
    //         case linelem::C: {
    //             auto e_c = (capacitor*)e;
    //             e_c->currents.push_back( 0 );
    //             E_i++;
    //             } break;
    //         default: break;
    //     }
    // }

}




void circuit::tran_fill_A(solver_t & A, Eigen::SparseMatrix<double> & A_copy, const size_t & n, const size_t & m) const {

    Eigen::MatrixXd A_dense, B, C, D, G;

    // Create matricies modeling storage elements using Norton companion model
    G = Eigen::MatrixXd::Zero(n, n);
    B = Eigen::MatrixXd::Zero(n, m);
    size_t E_i = 0;
    // for (size_t i = 0; i < (linelems.size()+itrelems.size()); i++) {
    //     const auto & e = (i<linelems.size()) ? linelems.at(i) : itrelems.at(i);
    for (const auto & e : linelems) {
        switch (e->ElemType) {
            case linelem::R: {
                auto e_r = (resistor*)e;
                if (e_r->Node1!=gnd)
                    G(e_r->Node1->i, e_r->Node1->i) += 1.0 / e_r->resistance;
                if (e_r->Node2!=gnd)
                    G(e_r->Node2->i, e_r->Node2->i) += 1.0 / e_r->resistance;
                if (e_r->Node1!=gnd && e_r->Node2!=gnd) {
                    G(e_r->Node1->i, e_r->Node2->i) -= 1.0 / e_r->resistance;
                    G(e_r->Node2->i, e_r->Node1->i) -= 1.0 / e_r->resistance;
                }
            } break;
            case linelem::C: {
                auto e_c = (capacitor*)e;
                if (e_c->Node1!=gnd) {
                    G(e_c->Node1->i, e_c->Node1->i) += e_c->conductance();
                }
                if (e_c->Node2!=gnd) {
                    G(e_c->Node2->i, e_c->Node2->i) += e_c->conductance();
                }
                if (e_c->Node1!=gnd && e_c->Node2!=gnd) {
                    G(e_c->Node1->i, e_c->Node2->i) -= e_c->conductance();
                    G(e_c->Node2->i, e_c->Node1->i) -= e_c->conductance();
                }
                } break;
            case linelem::L: {
                auto e_l = (inductor*)e;
                if (e_l->Node1!=gnd) {
                    G(e_l->Node1->i, e_l->Node1->i) += e_l->conductance();
                }
                if (e_l->Node2!=gnd) {
                    G(e_l->Node2->i, e_l->Node2->i) += e_l->conductance();
                }
                if (e_l->Node1!=gnd && e_l->Node2!=gnd) {
                    G(e_l->Node1->i, e_l->Node2->i) -= e_l->conductance();
                    G(e_l->Node2->i, e_l->Node1->i) -= e_l->conductance();
                }
                } break;
            case linelem::V: {
                auto e_v = (V_source*)e;
                if (e_v->Node1!=gnd)
                    B(e_v->Node1->i, E_i) = -1.0;
                if (e_v->Node2!=gnd)
                    B(e_v->Node2->i, E_i) = 1.0;
                E_i++;
                } break;
            default: break;
        }
    }

    C = B.transpose();
    D = Eigen::MatrixXd::Zero(m, m);
    A_dense = Eigen::MatrixXd(n+m, n+m);
    A_dense << G, B, C, D;

    A_copy = A_dense.sparseView();
    A.analyzePattern(A_copy);
    A.factorize(A_copy);

}



void circuit::tran_step(const solver_t & A, const Eigen::SparseMatrix<double> & A_copy, const size_t & n, const size_t & m) {

    // Generate z vector
    Eigen::Matrix <double, Eigen::Dynamic, 1> z, I, E;

    z = Eigen::MatrixXd::Zero(n+m, 1);
    I = Eigen::MatrixXd::Zero(n, 1);
    E = Eigen::MatrixXd::Zero(m, 1);

    size_t E_i = 0;
    for (auto e : linelems) {
        switch (e->ElemType) {
            case linelem::C: {
                auto e_c = (capacitor*)e;
                double I_eq = e_c->conductance() * e_c->voltage();
                if (e_c->Node1!=gnd)
                    I(e_c->Node1->i, 0) += I_eq;
                if (e_c->Node2!=gnd)
                    I(e_c->Node2->i, 0) -= I_eq;
                } break;
            case linelem::L: {
                auto e_l = (inductor*)e;
                double I_eq = e_l->current();
                if (e_l->Node1!=gnd)
                    I(e_l->Node1->i, 0) += I_eq;
                if (e_l->Node2!=gnd)
                    I(e_l->Node2->i, 0) -= I_eq;
                } break;
            case linelem::V: {
                auto e_v = (V_source*)e;
                E(E_i, 0) = -1.0 * e_v->voltage();
                E_i++;
                } break;
            case linelem::I: {
                auto e_i = (I_source*)e;
                if (e_i->Node1!=gnd)
                    I(e_i->Node1->i, 0) += e_i->current();
                if (e_i->Node2!=gnd)
                    I(e_i->Node2->i, 0) -= e_i->current();
                } break;
            default: break;
        }
    }

    Eigen::MatrixXd x;
    if (mosfets.size()==0) {
        z << I, E;
        x = A.solve(z);
        // cout << A_copy << endl << x << endl << z << endl << endl;
    } else {
        struct mosfet_nodes{double vd, vs, vg;};
        unordered_map<mosfet*,mosfet_nodes> k;
        for (auto m : mosfets) k[m] = {m->NodeD->voltage(), m->NodeS->voltage(), m->NodeG->voltage()};
        for (size_t nr_i = 0; nr_i < 100; nr_i++) {
            Eigen::SparseMatrix<double> NR_G(n+m, n+m);
            solver_t A_combined_solver;
            Eigen::Matrix <double, Eigen::Dynamic, 1> NR_I = Eigen::MatrixXd::Zero(n, 1);
            for (auto m : mosfets) {
                double V_GS = k[m].vg-k[m].vs;
                double V_DS = k[m].vd-k[m].vs;
                if (m->NodeS!=gnd) {
                    NR_I(m->NodeS->i, 0) += m->NR_I_eq(V_GS, V_DS);
                    NR_G.coeffRef(m->NodeS->i, m->NodeS->i) += m->NR_G_eq(V_GS, V_DS);
                }
                if (m->NodeD!=gnd) {
                    NR_I(m->NodeD->i, 0) -= m->NR_I_eq(V_GS, V_DS);
                    NR_G.coeffRef(m->NodeD->i, m->NodeD->i) += m->NR_G_eq(V_GS, V_DS);
                }
                if (m->NodeS!=gnd && m->NodeD!=gnd) {
                    NR_G.coeffRef(m->NodeD->i, m->NodeS->i) -= m->NR_G_eq(V_GS, V_DS);
                    NR_G.coeffRef(m->NodeS->i, m->NodeD->i) -= m->NR_G_eq(V_GS, V_DS);
                }
            }
            Eigen::SparseMatrix<double> A_combined = A_copy + NR_G;
            A_combined_solver.analyzePattern(A_combined);
            A_combined_solver.factorize(A_combined);
            z << (I+NR_I), E;
            while (A_combined_solver.info()) {
                cerr << "Program reached an unknown state. Exiting..." << endl;
                exit(1);
                // for (size_t i = 0; i < n+m; i++) {
                //     bool is_zero = true;
                //     for (size_t j = 0; j < n+m; j++) if (A_combined.coeffRef(i, j) != 0) is_zero = false;
                //     if (is_zero) {
                //         z(i,0) = 0;
                //         A_combined.coeffRef(i,0) = 1;
                //         A_combined.coeffRef(0,i) = 1;
                //         cout << "zero row found at " << i << endl;
                //     }
                // }
                // A_combined_solver.analyzePattern(A_combined);
                // A_combined_solver.factorize(A_combined);
            }

            x = A_combined_solver.solve(z);

            bool changed = false;
            for (auto m : mosfets) {
                double vs = k[m].vs;
                double vd = k[m].vd;
                double vg = k[m].vg;
                k[m].vs = ((m->NodeS==gnd) ? 0 : x(m->NodeS->i, 0));
                k[m].vd = ((m->NodeD==gnd) ? 0 : x(m->NodeD->i, 0));
                k[m].vg = ((m->NodeG==gnd) ? 0 : x(m->NodeG->i, 0));
                changed = changed || (abs(vs-k[m].vs)>tran_precision) || (abs(vd-k[m].vd)>tran_precision) || (abs(vg-k[m].vg)>tran_precision);
            }
            if (!changed) {
                // cout << "exited tran at nr_i=" << nr_i << endl;
                break;
            }
        }
    }


    // record voltages
    for (auto node_i : nodes) {
        // cout << "about to add " << x(node_i.second->i,0) << endl;
        node_i.second->voltages.push_back( x(node_i.second->i,0 ) );
    }

    // record currents
    E_i = 0;
    for (auto e : linelems) {
        switch (e->ElemType) {
            case linelem::C: {
                auto e_c = (capacitor*)e;
                double I = (e_c->conductance()*e_c->voltage(-2)) - (e_c->conductance()*e_c->voltage(-1));
                e_c->currents.push_back( I );
                } break;
            case linelem::L: {
                auto e_l = (inductor*)e;
                double I = e_l->current() - (e_l->conductance() * e_l->voltage(-1));
                e_l->currents.push_back( I );
                } break;
            case linelem::V: {
                auto e_v = (V_source*)e;
                e_v->currents.push_back( x(n+E_i, 0) );
                E_i++;
                } break;
            default: break;
        }
    }
    step_num++;
}


void circuit::tran() {

    // set initial voltages
    dc();

    solver_t A;
    Eigen::SparseMatrix<double> A_copy;

    // count number of each element
    size_t n, m;
    size_t num_capacitors = 0;
    size_t num_voltage_sources = 0;
    for (auto e : linelems) {
        switch (e->ElemType) {
            case linelem::C: num_capacitors++; break;
            case linelem::V: num_voltage_sources++; break;
            default: break;
        }
    }

    m = num_voltage_sources;
    n = nodes.size();

    circuit::tran_fill_A(A, A_copy, n, m);

    double time = 0;
    while (time < stop_time) {
        tran_step(A, A_copy, n, m);
        time = step_num * time_step;
    }
    stop_time = time;

}


// node getters
double circuit::node::voltage(const int & t) const {
    if (this==circuit::gnd)
        return 0;
    size_t i = (t<0) ? (voltages.size()+t) : (t);
    if (i >= voltages.size()) {
        cerr << "Warning: Tried to access voltage out of range." << endl;
        return NAN;
    }
    return voltages.at(i);
}

double circuit::mosfet::current(const size_t & t) const {
    double V_GS = NodeG->voltage(t) - NodeS->voltage(t);
    double V_DS = NodeD->voltage(t) - NodeS->voltage(t);
    return I_DS(V_GS, V_DS);
}


double circuit::mosfet::conductance(const double & V_GS, const double & V_DS) const {
    if (ElemType==nmos) {
        // cutoff
        if (V_GS <= V_T)
            return 0;

        // extra case for correct convergence
        if (V_DS < 0)
            return 1;

        double V_GST = V_GS-V_T;
        double k = MU * C_OX * W / L;
        double gm = 1 + LAMBDA * V_DS;

        // linear
        if (V_DS <= V_GST)
            return k*LAMBDA*(V_GST*V_DS-pow(V_DS,2)/2) + k*(-V_DS+V_GST)*gm;

        // saturation
        return 0.5 * k * pow(V_GST,2) * LAMBDA;
    } else if (ElemType==pmos) {
        // cutoff
        if (V_GS >= V_T)
            return 0;

        // extra case for correct convergence
        if (V_DS > 0)
            return 1;

        double V_GST = V_GS-V_T;
        double k = MU * C_OX * W / L;
        double gm = 1 - LAMBDA * V_DS;

        // linear
        if (V_DS >= V_GST)
            return k*LAMBDA*(V_GST*V_DS-pow(V_DS,2)/2) - k*(V_GST-V_DS)*gm;

        // saturation
        return 0.5 * k * pow(V_GST,2) * LAMBDA;
    } else {
        cerr << "Unknown MOSFET type: " << ElemType << endl;
        exit(1);
    }
}
double circuit::mosfet::I_DS(const double & V_GS, const double & V_DS) const {
    if (ElemType==nmos) {
        // cutoff
        if (V_GS <= V_T)
            return 0;

        // extra case for correct convergence
        if (V_DS < 0)
            return V_DS;

        double V_GST = V_GS-V_T;
        double k = MU * C_OX * W / L;
        double gm = 1 + (LAMBDA*V_DS);

        // linear
        if (V_DS <= V_GST) {
            return k * (V_DS*V_GST - 0.5*pow(V_DS,2)) * gm;
        }

        // saturation
        return 0.5 * k * pow(V_GST,2) * gm;
    } else if (ElemType==pmos) {
        // cutoff
        if (V_GS >= V_T)
            return 0;

        // extra case for correct convergence
        if (V_DS > 0)
            return V_DS;

        double V_GST = V_GS-V_T;
        double k = MU * C_OX * W / L;
        double gm = 1 - LAMBDA * V_DS;

        // linear
        if (V_DS >= V_GST)
            return -k * (V_DS*V_GST - 0.5*pow(V_DS,2)) * gm;

        // saturation
        return -0.5 * k * pow(V_GST,2) * gm;
    } else {
        cerr << "Unknown MOSFET type: " << ElemType << endl;
        exit(1);
    }
}

double circuit::mosfet::NR_G_eq(const double & V_GS,const double & V_DS) const {
    return conductance(V_GS, V_DS);
}

double circuit::mosfet::NR_I_eq(const double & V_GS,const double & V_DS) const {
    return I_DS(V_GS, V_DS) - (conductance(V_GS, V_DS)*V_DS);
}

// linelem const getters
double circuit::linelem::voltage(const int & t) const {
    return Node1->voltage(t) - Node2->voltage(t);
}
double circuit::resistor::current(const int & t) const {
    return voltage(t) / resistance;
}

// storage devices
double circuit::storage_device::current(const int & t) const {
    if (currents.size()==0) return NAN;
    return currents.at( (t<0)?(c.step_num+1+t):(t) );
}
double circuit::capacitor::conductance() const {
    return capacitance / c.time_step;
}
double circuit::inductor::conductance() const {
    return c.time_step / inductance;
}

// power sources
double circuit::V_source::current(const int & t) const {
    if (currents.size()==0) return NAN;
    return currents.at( (t<0)?(c.step_num+1+t):(t) );
}
double circuit::V_dc::voltage(const int & t) const {
    return voltage_value;
}
double pwl_value(double x, vector< pair<double,double> > v) {
    if (x < 0)
        return v.front().second;
    if (x > v.back().first)
        return v.back().second;
    for (size_t i = 1; i < v.size(); i++) {
        if (x < v.at(i).first)
            return v.at(i-1).second;
    }
    return v.at(0).second;
}
double circuit::V_pwl::voltage(const int & t) const {
    double time = (double)((t<0)?(c.step_num+1+t):(t)) * c.time_step;
    return pwl_value(time, voltages);
}
double circuit::I_dc::current(const int & t) const {
    return current_value;
}
double circuit::I_pwl::current(const int & t) const {
    double time = (double)((t<0)?(c.step_num+1+t):(t)) * c.time_step;
    return pwl_value(time, currents);
}


// print methods
void circuit::print() const {
    for (const auto & e : linelems)
        e->print();
    for (const auto & m : mosfets)
        m->print();
    for (const auto & n : nodes)
        n.second->print();
}

void circuit::node::print() const {
    cout << "n" << name << " id" << id << " i" << i << " V=" << voltage() << endl;
}

void circuit::linelem::print() const {
    cout << name << " n" << Node1->name << " n" << Node2->name << " V=" << voltage() << endl;
}
void circuit::resistor::print() const {
    cout << name << " n" << Node1->name << " n" << Node2->name << " R=" << resistance << " I=" << current() << " V=" << voltage() << endl;
}
void circuit::capacitor::print() const {
    cout << name << " n" << Node1->name << " n" << Node2->name << " C=" << capacitance << " V_i=" << initial_voltage << " I=" << current() << " V=" << voltage() << endl;
}
void circuit::inductor::print() const {
    cout << name << " n" << Node1->name << " n" << Node2->name << " L=" << inductance << " I_i=" << initial_current << " I=" << current() << " V=" << voltage() << endl;
}
void circuit::power_source::print() const {
    cout << name << " n" << Node1->name << " n" << Node2->name << " V=" << voltage() << endl;
}
void circuit::V_source::print() const {
    cout << name << " n" << Node1->name << " n" << Node2->name << " V=" << voltage() << " I=" << current() << endl;
}
void circuit::I_source::print() const {
    cout << name << " n" << Node1->name << " n" << Node2->name << " V=" << voltage() << " I=" << current() << endl;
}

void circuit::mosfet::print() const {
    cout << name << ((ElemType==nmos)?(" nmos"):(ElemType==pmos)?(" pmos"):(" ????")) << " D=n" << NodeD->name << " S=n" << NodeS->name << " G=n" << " I_DS=" << I_DS(NodeG->voltage()-NodeS->voltage(), NodeD->voltage()-NodeS->voltage()) << endl;
    // cout << name << ((ElemType==nmos)?(" nmos"):(ElemType==pmos)?(" pmos"):(" ????")) << " n" << Node1->name << " n" << Node2->name << " g_n" << NodeG->name << " W=" << W << " L=" << L << " V_T=" << V_T << " MU=" << MU << " C_OX=" << C_OX << " LAMBDA=" << LAMBDA << " C_J=" << C_J << endl;
}

void circuit::to_json(const std::string & filename) const {
    circuit_interface::export_circuit(this, filename);
}
