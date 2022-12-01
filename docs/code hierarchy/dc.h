
#include "analysis.h"
#include "circuit.h"
#include <Eigen/Eigen>
#include <unordered_map>

class dc : public analysis {
public:
    dc(circuit * c);

    void plotnv(matlab * m, int node_name) const;
    void printnv(int node_name) const;

    double voltage(circuit::node * n) const;
    double voltage(circuit::linelem * e) const;
    double current(circuit::linelem * e) const;

private:
    std::unordered_map< circuit::node*, double > node_voltage;
    std::unordered_map< circuit::capacitor*, double > capacitor_current;
    std::unordered_map< circuit::V_source*, double > V_source_current;

    struct mosfet_solver {
        static double conductance(circuit::mosfet * m, double V_GS, double V_DS);
        static double I_DS(circuit::mosfet * m, double V_GS, double V_DS);
        static double NR_G_eq(circuit::mosfet * m, double V_GS,const double V_DS);
        static double NR_I_eq(circuit::mosfet * m, double V_GS,const double V_DS);
    };

    void stamp_A(Eigen::SparseMatrix<double> A, circuit::linelem * e) const;
    void stamp_z(Eigen::SparseMatrix<double> z, circuit::linelem * e) const;
    void stamp_A_NR(Eigen::SparseMatrix<double> A_NR, circuit::mosfet * m, double V_GS, double V_DS) const;
    void stamp_z_NR(Eigen::SparseMatrix<double> z_NR, circuit::mosfet * m, double V_GS, double V_DS) const;
};
