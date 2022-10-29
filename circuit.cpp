
#include "circuit.h"
#include "circuit_loader.h"
#include <iostream>

using namespace std;

circuit::node GND_INSTANCE{.id=-1, .name=0, .i=-1, .voltage={0}};
circuit::node * const circuit::gnd = &GND_INSTANCE;

const double circuit::time_step = 0.0000001;

circuit::circuit() { }

circuit::circuit(const std::string & filename) {

    // initialize this->nodes and this->linelems
    circuit_loader::circuit_from_filename(this, filename);

    // init time
    time = 0;

    // count number of each element
    size_t num_capacitors = 0;
    size_t num_voltage_sources = 0;
    for (auto e : linelems) {
        switch (e->ElemType) {
            case linelem::C: num_capacitors++; break;
            case linelem::V: num_voltage_sources++; break;
            default: break;
        }
    }

    Eigen::MatrixXd B, C, D, E,    G,    I, x, z;
    size_t E_i;

    // set G array with resistors only
    n = nodes.size();
    G = Eigen::MatrixXd::Zero(n, n);
    for (auto e : linelems) {
        auto e_r = (resistor*)e;
        if (e->ElemType == linelem::R) {
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
        auto e_c = (capacitor*)e;
        auto e_l = (inductor*)e;
        auto e_i = (I_dc*)e;
        auto e_r = (resistor*)e;
        auto e_v = (V_dc*)e;
        switch (e->ElemType) {
            case linelem::V:
                if (e_v->Node1!=gnd)
                    B(e_v->Node1->i, E_i) = -1.0;
                if (e_v->Node2!=gnd)
                    B(e_v->Node2->i, E_i) = 1.0;
                E(E_i, 0) = -1.0 * e_v->voltage;
                E_i++;
                break;
            case linelem::C:
                if (e_c->Node1!=gnd)
                    B(e_c->Node1->i, E_i) = -1.0;
                if (e_c->Node2!=gnd)
                    B(e_c->Node2->i, E_i) = 1.0;
                E(E_i, 0) = -1.0 * e_c->initial_voltage;
                E_i++;
                break;
            case linelem::I:
                if (e_i->Node1!=gnd)
                    I(e_i->Node1->i, 0) += e_i->current;
                break;
            case linelem::L:
                if (e_l->Node1!=gnd)
                    I(e_l->Node1->i, 0) += e_l->initial_current;
                break;
            default: break;
        }
    }
    C = B.transpose();
    D = Eigen::MatrixXd::Zero(m, m);
    A = Eigen::MatrixXd(n+m, n+m);
    A << G, B, C, D;
    z = Eigen::MatrixXd(n+m, 1);
    z << I, E;
    x = A.completeOrthogonalDecomposition().solve(z);

    // record initial voltages
    for (auto node_i : nodes) {
        node_i.second->voltage.push_back( x(node_i.second->i,0 ) );
        cout << node_i.second->name << " voltage is " << x(node_i.second->i,0 ) << endl;
    }

    // record initial currents
    E_i = 0;
    for (auto e : linelems) {
        auto e_c = (capacitor*)e;
        auto e_v = (V_dc*)e;
        auto e_l = (inductor*)e;
        auto e_i = (I_dc*)e;
        switch (e->ElemType) {
            case linelem::I:
                e_i->currents.push_back( e_i->current );
                break;
            case linelem::L:
                e_l->currents.push_back( e_l->initial_current );
                break;
            case linelem::V:
                e_v->currents.push_back( x(n+E_i, 0) );
                E_i++;
                break;
            case linelem::C:
                e_c->currents.push_back( x(n+E_i, 0) );
                E_i++;
                break;
            default: break;
        }
    }

    // Find matricies modeling storage elements using Norton companion model
    m = num_voltage_sources;
    B = Eigen::MatrixXd::Zero(n, m);
    I = Eigen::MatrixXd::Zero(n, 1);
    E = Eigen::MatrixXd::Zero(m, 1);
    E_i = 0;
    for (auto e : linelems) {
        auto e_c = (capacitor*)e;
        auto e_l = (inductor*)e;
        auto e_i = (I_dc*)e;
        auto e_r = (resistor*)e;
        auto e_v = (V_dc*)e;
        switch (e->ElemType) {
            case linelem::V:
                if (e_v->Node1!=gnd)
                    B(e_v->Node1->i, E_i) = -1.0;
                if (e_v->Node2!=gnd)
                    B(e_v->Node2->i, E_i) = 1.0;
                E(E_i, 0) = -1.0 * e_v->voltage;
                E_i++;
                break;
            case linelem::C:
                if (e_c->Node1!=gnd) {
                    G(e_c->Node1->i, e_c->Node1->i) += 2.0 * e_c->capacitance / circuit::time_step;
                }
                if (e_c->Node2!=gnd) {
                    G(e_c->Node2->i, e_c->Node2->i) += 2.0 * e_c->capacitance / circuit::time_step;
                    double conductance = 2.0 * e_c->capacitance / circuit::time_step;
                    I(e_c->Node2->i, 0) -= e_c->initial_voltage * conductance; // ?
                }
                if (e_c->Node1!=gnd && e_c->Node2!=gnd) {
                    G(e_c->Node1->i, e_c->Node2->i) -= 2.0 * e_c->capacitance / circuit::time_step;
                    G(e_c->Node2->i, e_c->Node1->i) -= 2.0 * e_c->capacitance / circuit::time_step;
                }
                break;
            case linelem::I:
                if (e_i->Node1!=gnd)
                    I(e_i->Node1->i, 0) += e_i->currents.back();
                break;
            case linelem::L:
                if (e_l->Node1!=gnd) {
                    G(e_l->Node1->i, e_l->Node1->i) += circuit::time_step / (2.0 * e_l->inductance);
                }
                if (e_l->Node2!=gnd) {
                    G(e_l->Node2->i, e_l->Node2->i) += circuit::time_step / (2.0 * e_l->inductance);
                    double conductance = circuit::time_step / (2.0 * e_l->inductance);
                    double initial_voltage = e_l->Node1->voltage.back() - e_l->Node2->voltage.back();
                    I(e_l->Node2->i, 0) -= initial_voltage * conductance; // ?
                }
                if (e_l->Node1!=gnd && e_l->Node2!=gnd) {
                    G(e_l->Node1->i, e_l->Node2->i) -= circuit::time_step / (2.0 * e_l->inductance);
                    G(e_l->Node2->i, e_l->Node1->i) -= circuit::time_step / (2.0 * e_l->inductance);
                }
                break;
            default: break;
        }
    }
    C = B.transpose();
    D = Eigen::MatrixXd::Zero(m, m);
    A = Eigen::MatrixXd(n+m, n+m);
    A << G, B, C, D;
    z = Eigen::MatrixXd(n+m, 1);
    z << I, E;

    // verify
    x = A.completeOrthogonalDecomposition().solve(z);
    for (auto node_i : nodes)
        cout << node_i.second->name << " voltage is " << x(node_i.second->i,0 ) << endl;

}

void circuit::linelem::print() const {
    if (Node1 == 0 || Node2 == 0) {
        cerr << "Element Node not defined." << endl;
        exit(1);
    }
    switch (ElemType) {
        case C: cout << "C " << Node1->name << " " << Node2->name << endl; break;
        case I: cout << "I " << Node1->name << " " << Node2->name << endl; break;
        case L: cout << "L " << Node1->name << " " << Node2->name << endl; break;
        case R: cout << "R " << Node1->name << " " << Node2->name << endl; break;
        case V: cout << "V " << Node1->name << " " << Node2->name << endl; break;
        default: break;
    }
}

void circuit::node::print() const {
    cout
        << "Node - "
        << "id " << id
        << "; name " << name
        << endl;
}

void circuit::print() const {
    for (const auto & e : linelems)
        e->print();
    for (const auto & n : nodes)
        n.second->print();
    cout << A << endl;
}

circuit::~circuit() {
    for (auto e : linelems)
        if (e) delete e;
    for (auto n : nodes)
        if (n.second) delete n.second;
}
