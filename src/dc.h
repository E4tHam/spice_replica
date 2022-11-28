
#ifndef __DC_H
#define __DC_H

#include "analysis.h"

#include <unordered_map>

class dc : public analysis {
public:
    dc(const circuit * const c);

    double voltage(const circuit::node * const n) const {
        if (n==circuit::gnd) return 0;
        return node_voltage.at(n);
    }
    double voltage(const circuit::linelem * const e) const {
        switch (e->ElemType) {
            case circuit::linelem::V: return voltage((circuit::V_source*)e);
            case circuit::linelem::C: return voltage((circuit::capacitor*)e);
        }
        return voltage(e->Node1) - voltage(e->Node2);
    }
    double voltage(const circuit::capacitor * const c) const { return c->initial_voltage; }
    double voltage(const circuit::V_source * const v) const {
        switch (v->SOURCE_TYPE) {
            case circuit::power_source::DC: return voltage((circuit::V_dc*)v);
            case circuit::power_source::PWL: return voltage((circuit::V_pwl*)v);
            default: return NAN;
        }
    }
    double voltage(const circuit::V_dc * const v) const { return v->voltage_value; }
    double voltage(const circuit::V_pwl * const v) const { return v->voltages.front().second; }

    double current(const circuit::linelem * const e) const {
        switch (e->ElemType) {
            case circuit::linelem::R: return current((circuit::resistor*)e);
            case circuit::linelem::C: return current((circuit::capacitor*)e);
            case circuit::linelem::L: return current((circuit::inductor*)e);
            case circuit::linelem::V: return current((circuit::V_source*)e);
            case circuit::linelem::I: return current((circuit::I_source*)e);
            default: return NAN;
        }
    }
    double current(const circuit::resistor * const r) const { return voltage(r) / r->resistance; }
    double current(const circuit::capacitor * const c) const { return capacitor_current.at(c); }
    double current(const circuit::inductor * const l) const { return l->initial_current; }
    double current(const circuit::V_source * const v) const { return V_source_current.at(v); }
    double current(const circuit::I_source * const i) const {
        switch (i->SOURCE_TYPE) {
            case circuit::power_source::DC: return current((circuit::I_dc*)i);
            case circuit::power_source::PWL: return current((circuit::I_pwl*)i);
            default: return NAN;
        }
    }
    double current(const circuit::I_dc * const i) const { return i->current_value; }
    double current(const circuit::I_pwl * const i) const { return i->currents.front().second; }

    double current(const circuit::mosfet * const m) const {
        double V_GS = voltage(m->NodeG) - voltage(m->NodeS);
        double V_DS = voltage(m->NodeD) - voltage(m->NodeS);
        return mosfet_solver::I_DS(m, V_GS, V_DS);
    }

private:
    std::unordered_map< const circuit::node*, double > node_voltage;
    std::unordered_map< const circuit::capacitor*, double > capacitor_current;
    std::unordered_map< const circuit::V_source*, double > V_source_current;

    struct mosfet_solver {
        static double conductance(const circuit::mosfet * const m, const double & V_GS, const double & V_DS);
        static double I_DS(const circuit::mosfet * const m, const double & V_GS, const double & V_DS);
        static double NR_G_eq(const circuit::mosfet * const m, const double & V_GS,const double & V_DS);
        static double NR_I_eq(const circuit::mosfet * const m, const double & V_GS,const double & V_DS);
    };

    friend class tran;
};

#endif
