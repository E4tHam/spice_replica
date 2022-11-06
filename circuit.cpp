
#include "circuit.h"
#include "circuit_interface.h"
#include <iostream>
#include <math.h>

using namespace std;



// static variables
circuit::node GND_INSTANCE(-1, 0, -1);
circuit::node * const circuit::gnd = &GND_INSTANCE;





// circuit methods
circuit::circuit() { }

circuit::circuit(const std::string & filename) {

    // initialize this->nodes and this->linelems
    circuit_interface::circuit_from_filename(this, filename);

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
    Eigen::MatrixXd A_dense, B, C, D, G, x;
    Eigen::Matrix <double, Eigen::Dynamic, 1> I, E, z;
    size_t E_i;

    // nonlinear
    Eigen::MatrixXd NR_G;
    Eigen::Matrix <double, Eigen::Dynamic, 1> NR_I;

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
    z = Eigen::MatrixXd(n+m, 1);

    unordered_map<diode*,double> k;
    for (auto d : diodes) k[d] = 0.0;

    // https://www.desmos.com/calculator/z0wakcoerx nr
    // https://www.desmos.com/calculator/28kehu1od4 working math!
    for (size_t nr_i = 0; nr_i < 100; nr_i++) {
        NR_G = Eigen::MatrixXd::Zero(n, n);
        NR_I = Eigen::MatrixXd::Zero(n, 1);
        for (auto d : diodes) {
            if (d->Node1!=gnd) {
                NR_I(d->Node1->i, 0) -= d->NR_I_eq(k[d]);
                NR_G(d->Node1->i, d->Node1->i) += d->NR_G_eq(k[d]);
            }
            if (d->Node2!=gnd) {

                NR_I(d->Node2->i, 0) += d->NR_I_eq(k[d]);
                NR_G(d->Node2->i, d->Node2->i) += d->NR_G_eq(k[d]);
            }
            if (d->Node1!=gnd && d->Node2!=gnd) {
                NR_G(d->Node1->i, d->Node2->i) -= d->NR_G_eq(k[d]);
                NR_G(d->Node2->i, d->Node1->i) -= d->NR_G_eq(k[d]);
            }
        }
        A_dense << (G+NR_G), B, C, D;
        z << (I+NR_I), E;
        x = A_dense.completeOrthogonalDecomposition().solve(z);

        for (auto d : diodes) {
            double nr_voltage
                = ((d->Node1==gnd) ? 0 : x(d->Node1->i, 0))
                - ((d->Node2==gnd) ? 0 : x(d->Node2->i, 0));
            k[d] = nr_voltage;
        }
    }

    for (auto d : diodes) {
        d->voltages.push_back(k[d]);
    }

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

}




void circuit::tran_fill_A(solver_t & A, const size_t & n, const size_t & m) const {

    Eigen::MatrixXd A_dense, B, C, D, G;

    // Create matricies modeling storage elements using Norton companion model
    G = Eigen::MatrixXd::Zero(n, n);
    B = Eigen::MatrixXd::Zero(n, m);
    size_t E_i = 0;
    for (auto e : linelems) {
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
                double I_eq = e_c->voltage() * e_c->conductance();
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
                double I_eq = e_l->current();
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
                } break;
            default: break;
        }
    }

    C = B.transpose();
    D = Eigen::MatrixXd::Zero(m, m);
    A_dense = Eigen::MatrixXd(n+m, n+m);
    A_dense << G, B, C, D;

    A.analyzePattern(A_dense.sparseView());
    A.factorize(A_dense.sparseView());

}



void circuit::tran_step(const solver_t & A, const size_t & n, const size_t & m) {


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



    z << I, E;
    Eigen::MatrixXd x = A.solve(z);


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

    circuit::tran_fill_A(A, n, m);

    double time = 0;
    while (time < stop_time) {
        tran_step(A, n, m);
        time = step_num * time_step;
    }
    stop_time = time;

}


// node getters
double circuit::node::voltage(const int & t) const {
    if (this==circuit::gnd)
        return 0;
    else if (t<1)
        return voltages.at(voltages.size()+t);
    else
        return voltages.at(t);
}

double circuit::diode::NR_G_eq(const double & V) const {
    double I_0 = 0.0000000001;
    double V_T = 0.0259;
    return I_0 / V_T * exp(V/V_T);
}
double circuit::diode::NR_I_eq(const double & V) const {
    return I_D(V) - NR_G_eq(V)*V;
}
double circuit::diode::I_D(const double & V) const {
    double I_0 = 0.0000000001;
    double V_T = 0.0259;
    return I_0 * (exp(V/V_T)-1);
}

// linelem const getters
double circuit::linelem::voltage(const int & t) const {
    return Node1->voltage() - Node2->voltage();
}
double circuit::resistor::current(const int & t) const {
    return voltage(t) / resistance;
}

// storage devices
double circuit::storage_device::current(const int & t) const {
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
    for (const auto & d : diodes)
        d->print();
    for (const auto & n : nodes)
        n.second->print();
}

void circuit::node::print() const {
    cout << "n" << name << " V=" << voltage() << endl;
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
    cout << name << " n" << Node1->name << " n" << Node2->name << " I=" << current() << " I=" << current() << endl;
}

void circuit::diode::print() const {
    cout << "Diode n" << Node1->name << " n" << Node2->name << endl;
}

void circuit::to_json(const std::string & filename) const {
    circuit_interface::export_circuit(this, filename);
}
