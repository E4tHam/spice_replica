
#include "circuit_loader.h"
#include "json.hpp"
#include <iostream>
#include <fstream>

using namespace std;
using namespace nlohmann;

// set c->nodes and c->linelems
void circuit_loader::circuit_from_filename(circuit * c, const std::string & filename) {
    ifstream ifs(filename);
    auto j = json::parse( ifs );
    ifs.close();
    const auto LINELEM = j.at("LINELEM");
    const auto NODES = j.at("NODES");
    for (const auto & e_array : LINELEM) {

        circuit::linelem * e;

        // set ElemType
        auto ElemType = (circuit::linelem::ElemType_t)e_array.at(0);

        // set nodes (If node doesn't exist, add it to c->nodes)
        int node1_id = (int)e_array.at(2);
        if (    !(node1_id==circuit::gnd->id)
                && (c->nodes.find(node1_id)==c->nodes.end())
        ) {
            c->nodes[node1_id] = new circuit::node();
            c->nodes[node1_id]->id = node1_id;
            c->nodes[node1_id]->i = node1_id-1;
            try { // if only 1 element in NODES
                c->nodes[node1_id]->name = NODES.at(node1_id-1);
            } catch (nlohmann::json_abi_v3_11_2::detail::type_error e) {
                c->nodes[node1_id]->name = NODES;
            }
        }
        circuit::node *Node1 = (node1_id==circuit::gnd->id) ? circuit::gnd : c->nodes[node1_id];
        int node2_id = (int)e_array.at(3);
        if (    !(node2_id==circuit::gnd->id)
                && (c->nodes.find(node2_id)==c->nodes.end())
        ) {
            c->nodes[node2_id] = new circuit::node();
            c->nodes[node2_id]->id = node2_id;
            c->nodes[node2_id]->i = node2_id-1;
            try { // if only 1 element in NODES
                c->nodes[node2_id]->name = NODES.at(node2_id-1);
            } catch (nlohmann::json_abi_v3_11_2::detail::type_error e) {
                c->nodes[node2_id]->name = NODES;
            }
        }
        circuit::node *Node2 = (node2_id==circuit::gnd->id) ? circuit::gnd : c->nodes[node2_id];


        switch (ElemType) {
            case circuit::linelem::C:
                e = new circuit::capacitor(
                    ElemType,
                    Node1,
                    Node2,
                    (double)e_array.at(1),
                    (double)e_array.at(4)
                );
                break;
            case circuit::linelem::L:
                e = new circuit::inductor(
                    ElemType,
                    Node1,
                    Node2,
                    (double)e_array.at(1),
                    (double)e_array.at(4)
                );
                break;
            case circuit::linelem::V:
                switch ((circuit::power_source::V_TYPE_t)e_array.at(4)) {
                    case circuit::power_source::DC:
                        e = new circuit::V_dc(
                            ElemType,
                            Node1,
                            Node2,
                            (circuit::power_source::V_TYPE_t)e_array.at(4),
                            (double)e_array.at(1)
                        );
                        break;
                    // case circuit::power_source::PWL:
                    //     break;
                    default:
                        cerr << "Unknown V_TYPE: " << (circuit::power_source::V_TYPE_t)e_array.at(4) << endl;
                        exit(1);
                        break;
                }
                break;
            case circuit::linelem::I:
                switch ((circuit::power_source::V_TYPE_t)e_array.at(4)) {
                    case circuit::power_source::DC:
                        e = new circuit::I_dc(
                            ElemType,
                            Node1,
                            Node2,
                            (circuit::power_source::V_TYPE_t)e_array.at(4),
                            (double)e_array.at(1)
                        );
                        break;
                    // case circuit::power_source::PWL:
                    //     break;
                    default:
                        cerr << "Unknown V_TYPE: " << (circuit::power_source::V_TYPE_t)e_array.at(4) << endl;
                        exit(1);
                        break;
                }
                break;
            case circuit::linelem::R:
                e = new circuit::resistor(
                    ElemType,
                    Node1,
                    Node2,
                    (double)e_array.at(1)
                );
                break;
            default:
                cerr << "cl Unknown ElemType: " << ElemType << endl;
                exit(1);
                break;
        }

        c->linelems.push_back(e);
    }
}
