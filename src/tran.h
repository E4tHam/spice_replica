
#ifndef __TRAN_H
#define __TRAN_H

#include "helper.h"
#include "analysis.h"
#include "matlab.h"
#include <Eigen/Eigen>

#include <unordered_map>

class tran : public analysis {
public:
    tran(const circuit * const c, const double & time_step, const double & stop_time);

    void plotnv(matlab * const m, const int & node_name) const;
    void printnv(const int & node_name) const;

    double voltage(const circuit::node * const n, const int & t = -1) const {
        if (n==circuit::gnd) return 0;
        if (node_voltage.find(n) == node_voltage.end()) return NAN;
        return helper::bidirectional_access(node_voltage.at(n), t);
    }
    double voltage(const circuit::linelem * const e, const int & t = -1) const {
        if (e->ElemType==circuit::linelem::V) return voltage((circuit::V_source*)e, t);
        return voltage(e->Node1, t) - voltage(e->Node2, t);
    }
    double voltage(const circuit::V_source * const v, const int & t = -1) const {
        switch (v->SOURCE_TYPE) {
            case circuit::power_source::DC: return voltage((circuit::V_dc*)v, t);
            case circuit::power_source::PWL: return voltage((circuit::V_pwl*)v, t);
            default: return NAN;
        }
    }
    double voltage(const circuit::V_dc * const v, const int & t = -1) const { return v->voltage_value; }
    double voltage(const circuit::V_pwl * const v, const int & t = -1) const {
        double time = (t<0) ? (stop_time-(time_step*(t+1))) : (time_step*t);
        return helper::pwl_value(time, v->voltages);
    }
    double current(const circuit::linelem * const e, const int & t = -1) const {
        switch (e->ElemType) {
            case circuit::linelem::R: return current((circuit::resistor*)e, t);
            case circuit::linelem::C: return current((circuit::storage_device*)e, t);
            case circuit::linelem::L: return current((circuit::storage_device*)e, t);
            case circuit::linelem::V: return current((circuit::V_source*)e, t);
            case circuit::linelem::I: return current((circuit::I_source*)e, t);
            default: return NAN;
        }
    }
    double current(const circuit::resistor * const r, const int & t = -1) const { return voltage(r, t) / r->resistance; }
    double current(const circuit::storage_device * const s, const int & t = -1) const {
        return helper::bidirectional_access(storage_device_current.at(s), t);
    }
    double current(const circuit::V_source * const v, const int & t = -1) const {
        return helper::bidirectional_access(V_source_current.at(v), t);
    }
    double current(const circuit::I_source * const i, const int & t = -1) const {
        switch (i->SOURCE_TYPE) {
            case circuit::power_source::DC: return current((circuit::I_dc*)i, t);
            case circuit::power_source::PWL: return current((circuit::I_pwl*)i, t);
            default: return NAN;
        }
    }
    double current(const circuit::I_dc * const i, const int & t = -1) const { return i->current_value; }
    double current(const circuit::I_pwl * const i, const int & t = -1) const {
        return helper::pwl_value((time_step*t), i->currents);
    }

    double current(const circuit::mosfet * const m, const int & t = -1) const {
        double V_GS = voltage(m->NodeG, t) - voltage(m->NodeS, t);
        double V_DS = voltage(m->NodeD, t) - voltage(m->NodeS, t);
        return mosfet_solver::I_DS(m, V_GS, V_DS);
    }

private:
    std::unordered_map< const circuit::V_source*, size_t > V_source_i;

    std::unordered_map< const circuit::node*, std::vector<double> > node_voltage;
    std::unordered_map< const circuit::storage_device*, std::vector<double> > storage_device_current;
    std::unordered_map< const circuit::V_source*, std::vector<double> > V_source_current;

    struct storage_device_solver {
        static double conductance(const circuit::capacitor * const c, const double & time_step) {
            return c->capacitance / time_step;
        }
        static double conductance(const circuit::inductor * const l, const double & time_step) {
            return time_step / l->inductance;
        }
    };
    struct mosfet_solver {
        static double conductance(const circuit::mosfet * const m, const double & V_GS, const double & V_DS);
        static double I_DS(const circuit::mosfet * const m, const double & V_GS, const double & V_DS);
        static double NR_G_eq(const circuit::mosfet * const m, const double & V_GS, const double & V_DS);
        static double NR_I_eq(const circuit::mosfet * const m, const double & V_GS, const double & V_DS);
    };
    void stamp_A(Eigen::SparseMatrix<double> & A, const circuit::linelem * const e) const;
    void stamp_z(Eigen::SparseMatrix<double> & z, const circuit::linelem * const e) const;
    void extract_current(Eigen::SparseMatrix<double> & x, const circuit::linelem * const e);
    void stamp_A_NR(Eigen::SparseMatrix<double> & A_NR, const circuit::mosfet * const m, const double & V_GS, const double & V_DS) const;
    void stamp_z_NR(Eigen::SparseMatrix<double> & z_NR, const circuit::mosfet * const m, const double & V_GS, const double & V_DS) const;
    double time_step;
    double stop_time;

    friend class helper;
};

#endif
