
#include "helper.h"
#include <fstream>
#include <iostream>

using namespace std;

double helper::pwl_value(const double & x, std::vector< std::pair<double,double> > v) {
    if (x < 0)
        return v.front().second;
    if (x > v.back().first)
        return v.back().second;
    for (size_t i = 1; i < v.size(); i++) {
        if (x < v.at(i).first)
            return v.at(i-1).second;
    }
    return v.back().second;
}

double helper::bidirectional_access(const std::vector<double> & v, const int & bi) {
    size_t i = (bi<0) ? (bi + v.size()) : (bi);
    if (i >= v.size()) return NAN;
    return v.at(i);
}

void helper::export_circuit(const tran * const a, const std::string & filename) {
    using namespace nlohmann;
    json NODES = {};
    json PLOTNV = {};

    for (const auto & n : a->c->nodes) {
        NODES.push_back({
            {"name", n.second->name},
            {"voltages", a->node_voltage.at(n.second)}
        });
        PLOTNV.push_back(n.second->name);
    }
    json LINELEMS = {};
    for (auto e : a->c->linelems) {
        vector<double> voltages;
        size_t i = 0;
        for (double t = 0; t <= a->stop_time; t+=a->time_step) {
            voltages.push_back(a->voltage(e, i++));
        }
        LINELEMS.push_back({
            {"name", e->name},
            {"voltages", voltages}
        });
    }
    json NLNELEMS = {};
    // for (auto m : a->c->mosfets) {
    //     vector<double> currents;
    //     for (const ) {
    //         currents.push_back(m->current(i));
    //     }
    //     NLNELEMS.push_back({
    //         {"name", m->name},
    //         {"currents", currents}
    //     });
    // }
    json j = {
        {"time_step", a->time_step},
        {"stop_time", a->stop_time},
        {"NODES", NODES},
        // {"LINELEMS", LINELEMS},
        // {"NLNELEMS", NLNELEMS},
        {"PLOTNV", PLOTNV}
        // {"PLOTBV", c->PLOTBV},
        // {"PLOTBI", a->c->PLOTBI}
    };
    ofstream ofs(filename);
    ofs << j;
    ofs.close();
}
