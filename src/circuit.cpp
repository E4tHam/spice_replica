
#include "circuit.h"
#include <nlohmann/json.hpp>
#include <iostream>
#include <algorithm>
#include "analysis.h"

using namespace std;
using namespace nlohmann;

// static variables
circuit::node GND_INSTANCE(-1, 0, -1);
circuit::node * const circuit::gnd = &GND_INSTANCE;

// circuit methods
circuit::circuit() : run_pointer(0) { }

circuit::circuit(const std::string & js) : run_pointer(0) {
    auto j = json::parse(js);
    const auto LINELEM = j.at("LINELEM");
    const auto LINNAME = j.at("LINNAME");
    const auto NLNELEM = j.at("NLNELEM");
    const auto NLNNAME = j.at("NLNNAME");
    const auto NODES = j.at("NODES");
    const auto INFO = j.at("INFO");
    analysis_type = new analysis::TYPE_t(INFO.at(0));
    time_step = INFO.at(1);
    stop_time = INFO.at(2);

    vector<int> PLOTNV = j.at("PLOTNV").get< std::vector<int> >();
    // vector<int> PLOTBV = j.at("PLOTBV").get< std::vector<int> >();
    // vector<int> PLOTBI;
    // for (const auto & b : j.at("PLOTBI")) PLOTBI.push_back(b.at(0));

    size_t e_name_i = 0;
    for (const auto & e_array : NLNELEM) {
        int nodeD_id = (int)e_array.at(2);
        if ((nodeD_id != circuit::gnd->id) && (nodes.find(nodeD_id) == nodes.end())) {
            int nodeD_i = nodeD_id-1;
            int nodeD_name = NODES.at(nodeD_i);
            nodes[nodeD_id] = new circuit::node(nodeD_id, nodeD_name, nodeD_i);
        }
        circuit::node * NodeD = (nodeD_id==circuit::gnd->id) ? (circuit::gnd) : (nodes[nodeD_id]);
        int nodeG_id = (int)e_array.at(3);
        if ((nodeG_id != circuit::gnd->id) && (nodes.find(nodeG_id) == nodes.end())) {
            int nodeG_i = nodeG_id-1;
            int nodeG_name = NODES.at(nodeG_i);
            nodes[nodeG_id] = new circuit::node(nodeG_id, nodeG_name, nodeG_i);
        }
        circuit::node * NodeG = (nodeG_id==circuit::gnd->id) ? (circuit::gnd) : (nodes[nodeG_id]);
        int nodeS_id = (int)e_array.at(4);
        if ((nodeS_id != circuit::gnd->id) && (nodes.find(nodeS_id) == nodes.end())) {
            int nodeS_i = nodeS_id-1;
            int nodeS_name = NODES.at(nodeS_i);
            nodes[nodeS_id] = new circuit::node(nodeS_id, nodeS_name, nodeS_i);
        }
        circuit::node * NodeS = (nodeS_id==circuit::gnd->id) ? (circuit::gnd) : (nodes[nodeS_id]);

        // set name
        string name = "";
        for (auto x : NLNNAME.at(e_name_i)) name += (char)(int)x;
        e_name_i++;

        mosfets.push_back(new circuit::mosfet(
            *this,
            (circuit::mosfet::ElemType_t)e_array.at(1),
            NodeD,
            NodeS,
            NodeG,
            name,
            (double)e_array.at(5),  // W
            (double)e_array.at(6),  // L
            (double)e_array.at(7),  // V_T
            (double)e_array.at(8),  // MU
            (double)e_array.at(9),  // C_OX
            (double)e_array.at(10), // LAMBDA
            (double)e_array.at(11)  // C_J
        ));

    }


    e_name_i = 0;
    for (const auto & e_array : LINELEM) {

        circuit::linelem * e;

        // set ElemType
        auto ElemType = (circuit::linelem::ElemType_t)e_array.at(0);

        // set name
        string name = "";
        for (auto x : LINNAME.at(e_name_i)) name += (char)(int)x;
        e_name_i++;

        // set nodes (If node doesn't exist, add it to nodes)
        int node1_id = (int)e_array.at(2);
        if (    !(node1_id==circuit::gnd->id)
                && (nodes.find(node1_id)==nodes.end())
        ) {
            int name;
            name = NODES.at(node1_id-1);
            nodes[node1_id] = new circuit::node(
                node1_id,
                name,
                node1_id-1
            );
        }
        circuit::node *Node1 = (node1_id==circuit::gnd->id) ? circuit::gnd : nodes[node1_id];
        int node2_id = (int)e_array.at(3);
        if (    !(node2_id==circuit::gnd->id)
                && (nodes.find(node2_id)==nodes.end())
        ) {
            int name;
            name = NODES.at(node2_id-1);
            nodes[node2_id] = new circuit::node(
                node2_id,
                name,
                node2_id-1
            );
        }
        circuit::node *Node2 = (node2_id==circuit::gnd->id) ? circuit::gnd : nodes[node2_id];

        // Parse e_array into new element
        switch (ElemType) {
            case circuit::linelem::R:
                e = new circuit::resistor(
                    *this,
                    ElemType,
                    name,
                    Node1,
                    Node2,
                    (double)e_array.at(1)
                );
                break;
            case circuit::linelem::C: {
                double initial_voltage = 0.0;
                try {initial_voltage = (double)e_array.at(4);}
                catch (nlohmann::json_abi_v3_11_2::detail::type_error e) { }
                e = new circuit::capacitor(
                    *this,
                    ElemType,
                    name,
                    Node1,
                    Node2,
                    (double)e_array.at(1),
                    initial_voltage
                );
                } break;
            case circuit::linelem::L: {
                double initial_current = 0.0;
                try {initial_current = (double)e_array.at(4);}
                catch (nlohmann::json_abi_v3_11_2::detail::type_error e) { }
                e = new circuit::inductor(
                    *this,
                    ElemType,
                    name,
                    Node1,
                    Node2,
                    (double)e_array.at(1),
                    initial_current
                );
                } break;
            case circuit::linelem::V:
                switch ((circuit::power_source::TYPE_t)e_array.at(4)) {
                    case circuit::power_source::DC: {
                        e = new circuit::V_dc(
                            *this,
                            ElemType,
                            name,
                            Node1,
                            Node2,
                            (circuit::power_source::TYPE_t)e_array.at(4),
                            (double)e_array.at(1)
                        );
                        } break;
                    case circuit::power_source::PWL: {
                        circuit::V_pwl::voltages_t voltages = {{0, e_array.at(1)}};
                        for (int i = 6; i < e_array.size(); i+=2)
                            voltages.emplace_back(e_array.at(i), e_array.at(i+1));
                        e = new circuit::V_pwl(
                            *this,
                            ElemType,
                            name,
                            Node1,
                            Node2,
                            (circuit::power_source::TYPE_t)e_array.at(4),
                            voltages
                        );
                        } break;
                    default:
                        cerr << "Unknown SOURCE_TYPE: " << (circuit::power_source::TYPE_t)e_array.at(4) << endl;
                        throw circuit::power_source::UnsupportedPowerSourceType();
                        break;
                }
                break;
            case circuit::linelem::I:
                switch ((circuit::power_source::TYPE_t)e_array.at(4)) {
                    case circuit::power_source::DC:
                        e = new circuit::I_dc(
                            *this,
                            ElemType,
                            name,
                            Node1,
                            Node2,
                            (circuit::power_source::TYPE_t)e_array.at(4),
                            (double)e_array.at(1)
                        );
                        break;
                    case circuit::power_source::PWL: {
                        circuit::I_pwl::currents_t currents = {{0, e_array.at(1)}};
                        for (int i = 6; i < e_array.size(); i+=2)
                            currents.emplace_back(e_array.at(i), e_array.at(i+1));
                        e = new circuit::I_pwl(
                            *this,
                            ElemType,
                            name,
                            Node1,
                            Node2,
                            (circuit::power_source::TYPE_t)e_array.at(4),
                            currents
                        );
                        } break;
                    default:
                        cerr << "Unknown SOURCE_TYPE: " << (circuit::power_source::TYPE_t)e_array.at(4) << endl;
                        throw circuit::power_source::UnsupportedPowerSourceType();
                        break;
                }
                break;
            default:
                throw circuit::linelem::UnsupportedLinearElementType();
                break;
        }

        linelems.push_back(e);
    }
    for (auto n : PLOTNV) {
        this->PLOTNV.push_back(nodes.at(n)->name);
    }
    std::reverse(this->PLOTNV.begin(), this->PLOTNV.end());
    // for (auto n : PLOTBV) PLOTBV.push_back(nodes.at(n)->name);
    // for (auto n : PLOTBI) PLOTBI.push_back(nodes.at(n)->name);
}

circuit::~circuit() {
    for (auto e : linelems)
        if (e) delete e;
    for (auto n : nodes)
        if (n.second) delete n.second;
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
    cout << "n" << name << " id" << id << " i" << i << endl;
}

void circuit::linelem::print() const {
    cout << name << " n" << Node1->name << " n" << Node2->name << endl;
}
void circuit::resistor::print() const {
    cout << name << " n" << Node1->name << " n" << Node2->name << endl;
}
void circuit::capacitor::print() const {
    cout << name << " n" << Node1->name << " n" << Node2->name << endl;
}
void circuit::inductor::print() const {
    cout << name << " n" << Node1->name << " n" << Node2->name << endl;
}
void circuit::power_source::print() const {
    cout << name << " n" << Node1->name << " n" << Node2->name << endl;
}
void circuit::V_source::print() const {
    cout << name << " n" << Node1->name << " n" << Node2->name << endl;
}
void circuit::I_source::print() const {
    cout << name << " n" << Node1->name << " n" << Node2->name << endl;
}
void circuit::mosfet::print() const {
    cout << name << ((ElemType==nmos)?(" nmos"):(ElemType==pmos)?(" pmos"):(" ????")) << " ND=" << NodeD->name << " NG=" << NodeG->name << " NS=" << NodeS->name << " W=" << W << " L=" << L << " V_T=" << V_T << " MU=" << MU << " C_OX=" << C_OX << " LAMBDA=" << LAMBDA << " C_J=" << C_J << endl;
}

analysis * circuit::run(matlab * const m) const {
    if (run_pointer!=0)
        delete run_pointer;
    switch (*((analysis::TYPE_t*)analysis_type)) {
        case analysis::DC: {
            run_pointer = new dc(this);
            return run_pointer;
            } break;
        case analysis::TRAN_FE:
        case analysis::TRAN_BE:
        case analysis::TRAN_TR: {
            run_pointer = new tran(this, time_step, stop_time);
            for (const auto & n : PLOTNV) {
                run_pointer->plotnv(m, n);
            }
            return run_pointer;
            } break;
        default: throw analysis::UnsupportedAnalysisType();
    }
}
