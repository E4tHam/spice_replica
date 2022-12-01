
#include "helper.h"
#include "analysis.h"
#include "matlab.h"
#include <Eigen/Eigen>
#include <unordered_map>

class tran : public analysis {
public:
    tran(circuit * c, double time_step, double stop_time);

    void plotnv(matlab * m, int node_name) const;
    void printnv(int node_name) const;

    double voltage(circuit::node * n, int t);
    double voltage(circuit::linelem * e, int t) const;
    double current(circuit::linelem * e, int t) const;

private:
    std::unordered_map< circuit::node*, std::vector<double> > node_voltage;
    std::unordered_map< circuit::storage_device*, std::vector<double> > storage_device_current;
    std::unordered_map< circuit::V_source*, std::vector<double> > V_source_current;

    struct storage_device_solver {
        static double conductance(circuit::capacitor * c, double time_step);
        static double conductance(circuit::inductor * l, double time_step);
    };
    struct mosfet_solver {
        static double conductance(circuit::mosfet * m, double V_GS, double V_DS);
        static double I_DS(circuit::mosfet * m, double V_GS, double V_DS);
        static double NR_G_eq(circuit::mosfet * m, double V_GS, double V_DS);
        static double NR_I_eq(circuit::mosfet * m, double V_GS, double V_DS);
    };
    void stamp_A(Eigen::SparseMatrix<double> A, circuit::linelem * e) const;
    void stamp_z(Eigen::SparseMatrix<double> z, circuit::linelem * e) const;
    void extract_current(Eigen::SparseMatrix<double> x, circuit::linelem * e);
    void stamp_A_NR(Eigen::SparseMatrix<double> A_NR, circuit::mosfet * m, double V_GS, double V_DS) const;
    void stamp_z_NR(Eigen::SparseMatrix<double> z_NR, circuit::mosfet * m, double V_GS, double V_DS) const;
    double time_step;
    double stop_time;
};
