
#ifndef __CIRCUIT_H
#define __CIRCUIT_H

#include <string>
#include <vector>
#include <unordered_map>
#include <Eigen/Dense>

class circuit {
public:
    // structs
    struct node {
        int id, name, i;
        std::vector<double> voltages;
        void print() const;
    };
    struct linelem {
        enum ElemType_t {
            R = 82,
            C = 67,
            L = 76,
            V = 86,
            I = 73
        };
        ElemType_t ElemType;
        node *Node1, *Node2;
        linelem(ElemType_t ElemType, node *Node1, node *Node2)
            : ElemType(ElemType), Node1(Node1), Node2(Node2) { }
        virtual void print() const;
    };
    struct resistor : linelem {
        double resistance;
        resistor(ElemType_t ElemType, node *Node1, node *Node2, double resistance)
            : linelem(ElemType, Node1, Node2), resistance(resistance) { }
        void print() const;
    };
    struct capacitor : linelem {
        double capacitance, initial_voltage;
        std::vector<double> currents;
        capacitor(ElemType_t ElemType, node *Node1, node *Node2, double capacitance, double initial_voltage)
            : linelem(ElemType, Node1, Node2), capacitance(capacitance), initial_voltage(initial_voltage) { }
        void print() const;
    };
    struct inductor : linelem {
        double inductance, initial_current;
        std::vector<double> currents;
        inductor(ElemType_t ElemType, node *Node1, node *Node2, double inductance, double initial_current)
            : linelem(ElemType, Node1, Node2), inductance(inductance), initial_current(initial_current) { }
        void print() const;
    };
    struct power_source : linelem {
        enum V_TYPE_t {
            DC = 0,
            AC = 1,
            PWL = 2
        };
        V_TYPE_t V_TYPE;
        std::vector<double> currents;
        power_source(ElemType_t ElemType, node *Node1, node *Node2, V_TYPE_t V_TYPE)
            : linelem(ElemType, Node1, Node2), V_TYPE(V_TYPE) { }
        void print() const;
    };
    struct V_dc : power_source {
        double voltage;
        V_dc(ElemType_t ElemType, node *Node1, node *Node2, V_TYPE_t V_TYPE, double voltage)
            : power_source(ElemType, Node1, Node2, V_TYPE), voltage(voltage) { }
        void print() const;
    };
    struct I_dc : power_source {
        double current;
        I_dc(ElemType_t ElemType, node *Node1, node *Node2, V_TYPE_t V_TYPE, double current)
            : power_source(ElemType, Node1, Node2, V_TYPE), current(current) { }
        void print() const;
    };
    // struct V_pwl : power_source { };
    // struct I_pwl : power_source { };


    // static variables
    static node * const gnd;
    static const double time_step;



    // functions
    circuit();
    circuit(const std::string & filename);
    ~circuit();

    void print() const;

private:
    // time
    double time;

    // elements and nodes
    std::vector<linelem*> linelems;
    std::unordered_map<int,node*> nodes;

    // matrix model
    Eigen::MatrixXd A;
    Eigen::Matrix <double, Eigen::Dynamic, 1> z;
    size_t n, m;

    // friends
    friend class circuit_loader;
};

#endif
