
#ifndef __CIRCUIT_H
#define __CIRCUIT_H

#include <string>
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>

class circuit {
public:
    struct node {
        int id, name, i;
        std::vector<double> voltage;
        void print() const;
    };
    static node * const gnd;
    static const double time_step;
    struct linelem {
        void print() const;
        enum ElemType_t {
            C = 67,
            I = 73,
            L = 76,
            R = 82,
            V = 86
        };
        ElemType_t ElemType;
        node *Node1, *Node2;
        linelem(ElemType_t ElemType, node *Node1, node *Node2)
            : ElemType(ElemType), Node1(Node1), Node2(Node2) { }
    };

    struct capacitor : linelem {
        double capacitance, initial_voltage;
        std::vector<double> currents;
        capacitor(ElemType_t ElemType, node *Node1, node *Node2, double capacitance, double initial_voltage)
            : linelem(ElemType, Node1, Node2), capacitance(capacitance), initial_voltage(initial_voltage) { }
    };
    struct inductor : linelem {
        double inductance, initial_current;
        std::vector<double> currents;
        inductor(ElemType_t ElemType, node *Node1, node *Node2, double inductance, double initial_current)
            : linelem(ElemType, Node1, Node2), inductance(inductance), initial_current(initial_current) { }
    };

    struct power_source : linelem {
        enum V_TYPE_t {
            DC = 0,
            AC = 1,
            PWL = 2
        };
        V_TYPE_t V_TYPE;
        power_source(ElemType_t ElemType, node *Node1, node *Node2, V_TYPE_t V_TYPE)
            : linelem(ElemType, Node1, Node2), V_TYPE(V_TYPE) { }
    };
    struct I_dc : power_source {
        double current;
        std::vector<double> currents;
        I_dc(ElemType_t ElemType, node *Node1, node *Node2, V_TYPE_t V_TYPE, double current)
            : power_source(ElemType, Node1, Node2, V_TYPE), current(current) { }
    };
    struct V_dc : power_source {
        double voltage;
        std::vector<double> currents;
        V_dc(ElemType_t ElemType, node *Node1, node *Node2, V_TYPE_t V_TYPE, double voltage)
            : power_source(ElemType, Node1, Node2, V_TYPE), voltage(voltage) { }
    };
    // struct I_pwl : power_source { };
    // struct V_pwl : power_source { };

    struct resistor : linelem {
        double resistance;
        resistor(ElemType_t ElemType, node *Node1, node *Node2, double resistance)
            : linelem(ElemType, Node1, Node2), resistance(resistance) { }
    };

    circuit();
    circuit(const std::string & filename);
    void print() const;
    ~circuit();
private:
    double time;
    std::vector<linelem*> linelems;
    std::unordered_map<int,node*> nodes;
    friend class circuit_loader;
    Eigen::MatrixXd A;
    Eigen::Matrix <double, Eigen::Dynamic, 1> z;
    size_t n, m;
};

#endif
