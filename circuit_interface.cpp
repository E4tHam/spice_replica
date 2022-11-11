
#include "circuit_interface.h"
#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>

using namespace std;
using namespace nlohmann;

circuit::node * circuit_interface::add_node(circuit * c, int node_id, int name) {
        // set nodes (If node doesn't exist, add it to c->nodes)
        if ( c->nodes.find(node_id) == c->nodes.end() ) {
            c->nodes[node_id] = new circuit::node(
                node_id,
                name,
                node_id-1
            );
        }
        return (node_id < 0) ? circuit::gnd : c->nodes[node_id];
}

// set c->nodes and c->linelems
void circuit_interface::circuit_from_filename(circuit * c, const std::string & filename) {
    ifstream ifs(filename);
    auto j = json::parse( ifs );
    ifs.close();
    const auto LINELEM = j.at("LINELEM");
    const auto LINNAME = j.at("LINNAME");
    const auto NLNELEM = j.at("NLNELEM");
    const auto NLNNAME = j.at("NLNNAME");
    const auto NODES = j.at("NODES");
    const auto INFO = j.at("INFO");
    if (INFO.at(0) != 3)
        throw TranNotFound();
    c->time_step = INFO.at(1);
    c->stop_time = INFO.at(2);

    vector<int> PLOTNV;
    // vector<int> PLOTBV;
    vector<int> PLOTBI;
    try {PLOTNV = j.at("PLOTNV").get< std::vector<int> >();}
    catch (nlohmann::json_abi_v3_11_2::detail::type_error e) { PLOTNV = {j.at("PLOTNV")}; }
    // try {PLOTBV = j.at("PLOTBV").get< std::vector<int> >();}
    // catch (nlohmann::json_abi_v3_11_2::detail::type_error e) { PLOTBV = {j.at("PLOTBV")}; }
    try { for (const auto & b : j.at("PLOTBI")) PLOTBI.push_back(b.at(0)); }
    catch (nlohmann::json_abi_v3_11_2::detail::type_error e) { PLOTBI = {j.at("PLOTBI").at(0)}; }

    size_t e_name_i = 0;
    try {
        for (const auto & e_array : NLNELEM) {
            circuit::node * NodeD = circuit::gnd;
            int nodeD_id = (int)e_array.at(2);
            if (nodeD_id != circuit::gnd->id) {
                int nodeD_name;
                try {nodeD_name = NODES.at(nodeD_id-1);}
                catch (nlohmann::json_abi_v3_11_2::detail::type_error e) { nodeD_name = NODES; }
                NodeD = add_node(c, nodeD_id, nodeD_name);
            }
            circuit::node * NodeS = circuit::gnd;
            int nodeS_id = (int)e_array.at(4);
            if (nodeS_id != circuit::gnd->id) {
                int nodeS_name;
                try {nodeS_name = NODES.at(nodeS_id-1);}
                catch (nlohmann::json_abi_v3_11_2::detail::type_error e) { nodeS_name = NODES; }
                NodeS = add_node(c, nodeS_id, nodeS_name);
            }
            circuit::node * NodeG = circuit::gnd;
            int nodeG_id = (int)e_array.at(3);
            if (nodeG_id != circuit::gnd->id) {
                int nodeG_name;
                try {nodeG_name = NODES.at(nodeG_id-1);}
                catch (nlohmann::json_abi_v3_11_2::detail::type_error e) { nodeG_name = NODES; }
                NodeG = add_node(c, nodeG_id, nodeG_name);
            }

            // set name
            string name = "";
            for (auto x : NLNNAME.at(e_name_i)) name += (char)(int)x;
            e_name_i++;

            c->mosfets.push_back(new circuit::mosfet(
                *c,
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

            c->itrelems.push_back(new circuit::capacitor(
                *c,
                circuit::linelem::C,
                (name+"__C_GS"),
                NodeG,
                NodeS,
                (0.5 * (double)e_array.at(9) * (double)e_array.at(5) * (double)e_array.at(6)),
                NAN
            ));
            c->itrelems.push_back(new circuit::capacitor(
                *c,
                circuit::linelem::C,
                (name+"__C_GD"),
                NodeG,
                NodeD,
                (0.5 * (double)e_array.at(9) * (double)e_array.at(5) * (double)e_array.at(6)),
                NAN
            ));
            c->itrelems.push_back(new circuit::capacitor(
                *c,
                circuit::linelem::C,
                (name+"__C_D0"),
                NodeD,
                circuit::gnd,
                (double)e_array.at(11),
                NAN
            ));
            c->itrelems.push_back(new circuit::capacitor(
                *c,
                circuit::linelem::C,
                (name+"__C_S0"),
                NodeS,
                circuit::gnd,
                (double)e_array.at(11),
                NAN
            ));

        }
    } catch (nlohmann::json_abi_v3_11_2::detail::type_error e) {
        circuit::node * NodeD = circuit::gnd;
        int nodeD_id = (int)NLNELEM.at(2);
        if (nodeD_id != circuit::gnd->id) {
            int nodeD_name;
            try {nodeD_name = NODES.at(nodeD_id-1);}
            catch (nlohmann::json_abi_v3_11_2::detail::type_error e) { nodeD_name = NODES; }
            NodeD = add_node(c, nodeD_id, nodeD_name);
        }
        circuit::node * NodeS = circuit::gnd;
        int nodeS_id = (int)NLNELEM.at(4);
        if (nodeS_id != circuit::gnd->id) {
            int nodeS_name;
            try {nodeS_name = NODES.at(nodeS_id-1);}
            catch (nlohmann::json_abi_v3_11_2::detail::type_error e) { nodeS_name = NODES; }
            NodeS = add_node(c, nodeS_id, nodeS_name);
        }
        circuit::node * NodeG = circuit::gnd;
        int nodeG_id = (int)NLNELEM.at(3);
        if (nodeG_id != circuit::gnd->id) {
            int nodeG_name;
            try {nodeG_name = NODES.at(nodeG_id-1);}
            catch (nlohmann::json_abi_v3_11_2::detail::type_error e) { nodeG_name = NODES; }
            NodeG = add_node(c, nodeG_id, nodeG_name);
        }

        // set name
        string name = "";
        for (auto x : NLNNAME.at(e_name_i)) name += (char)(int)x;
        e_name_i++;

        c->mosfets.push_back(new circuit::mosfet(
            *c,
            (circuit::mosfet::ElemType_t)NLNELEM.at(1),
            NodeD,
            NodeS,
            NodeG,
            name,
            (double)NLNELEM.at(5),
            (double)NLNELEM.at(6),
            (double)NLNELEM.at(7),
            (double)NLNELEM.at(8),
            (double)NLNELEM.at(9),
            (double)NLNELEM.at(10),
            (double)NLNELEM.at(11)
        ));

        c->itrelems.push_back(new circuit::capacitor(
            *c,
            circuit::linelem::C,
            (name+"__C_GS"),
            NodeG,
            NodeS,
            (0.5 * (double)NLNELEM.at(9) * (double)NLNELEM.at(5) * (double)NLNELEM.at(6)),
            NAN
        ));
        c->itrelems.push_back(new circuit::capacitor(
            *c,
            circuit::linelem::C,
            (name+"__C_GD"),
            NodeG,
            NodeD,
            (0.5 * (double)NLNELEM.at(9) * (double)NLNELEM.at(5) * (double)NLNELEM.at(6)),
            NAN
        ));
        c->itrelems.push_back(new circuit::capacitor(
            *c,
            circuit::linelem::C,
            (name+"__C_D0"),
            NodeD,
            circuit::gnd,
            (double)NLNELEM.at(11),
            NAN
        ));
        c->itrelems.push_back(new circuit::capacitor(
            *c,
            circuit::linelem::C,
            (name+"__C_S0"),
            NodeS,
            circuit::gnd,
            (double)NLNELEM.at(11),
            NAN
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

        // set nodes (If node doesn't exist, add it to c->nodes)
        int node1_id = (int)e_array.at(2);
        if (    !(node1_id==circuit::gnd->id)
                && (c->nodes.find(node1_id)==c->nodes.end())
        ) {
            int name;
            try {name = NODES.at(node1_id-1);}
            catch (nlohmann::json_abi_v3_11_2::detail::type_error e) { name = NODES; }
            c->nodes[node1_id] = new circuit::node(
                node1_id,
                name,
                node1_id-1
            );
        }
        circuit::node *Node1 = (node1_id==circuit::gnd->id) ? circuit::gnd : c->nodes[node1_id];
        int node2_id = (int)e_array.at(3);
        if (    !(node2_id==circuit::gnd->id)
                && (c->nodes.find(node2_id)==c->nodes.end())
        ) {
            int name;
            try {name = NODES.at(node2_id-1);}
            catch (nlohmann::json_abi_v3_11_2::detail::type_error e) { name = NODES; }
            c->nodes[node2_id] = new circuit::node(
                node2_id,
                name,
                node2_id-1
            );
        }
        circuit::node *Node2 = (node2_id==circuit::gnd->id) ? circuit::gnd : c->nodes[node2_id];

        // Parse e_array into new element
        switch (ElemType) {
            case circuit::linelem::R:
                e = new circuit::resistor(
                    *c,
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
                    *c,
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
                    *c,
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
                            *c,
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
                            *c,
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
                        exit(1);
                        break;
                }
                break;
            case circuit::linelem::I:
                switch ((circuit::power_source::TYPE_t)e_array.at(4)) {
                    case circuit::power_source::DC:
                        e = new circuit::I_dc(
                            *c,
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
                            *c,
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
                        exit(1);
                        break;
                }
                break;
            default:
                cerr << "Unknown ElemType: " << ElemType << endl;
                exit(1);
                break;
        }

        c->linelems.push_back(e);
    }

    for (auto n : PLOTNV) c->PLOTNV.push_back(c->nodes.at(n)->name);
    // for (auto n : PLOTBV) c->PLOTBV.push_back(c->nodes.at(n)->name);
    for (auto n : PLOTBI) c->PLOTBI.push_back(c->nodes.at(n)->name);
}





void circuit_interface::export_circuit(const circuit * const c, const std::string & filename) {
    json NODES = {};
    for (auto n : c->nodes) {
        NODES.push_back({
            {"name", n.second->name},
            {"voltages", n.second->voltages}
        });
    }
    json LINELEMS = {};
    for (auto e : c->linelems) {
        vector<double> voltages;
        for (size_t i = 0; i <= c->step_num; i++) {
            voltages.push_back(e->voltage(i));
        }
        LINELEMS.push_back({
            {"name", e->name},
            {"voltages", voltages}
        });
    }
    json NLNELEMS = {};
    for (auto m : c->mosfets) {
        vector<double> currents;
        for (size_t i = 0; i <= c->step_num; i++) {
            currents.push_back(m->current(i));
        }
        NLNELEMS.push_back({
            {"name", m->name},
            {"currents", currents}
        });
    }
    json j = {
        {"time_step", c->time_step},
        {"stop_time", c->stop_time},
        {"NODES", NODES},
        {"LINELEMS", LINELEMS},
        {"NLNELEMS", NLNELEMS},
        {"PLOTNV", c->PLOTNV},
        // {"PLOTBV", c->PLOTBV},
        {"PLOTBI", c->PLOTBI}
    };
    ofstream ofs(filename);
    ofs << j;
    ofs.close();
}
