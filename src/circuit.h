
#ifndef __CIRCUIT_H
#define __CIRCUIT_H

#include <string>
#include <vector>
#include <unordered_map>
#include <Eigen/Eigen>
#include <utility>
#include "matlab.h"
#include <memory>

class analysis;

class circuit {
public:
    // structs
    struct node {
        int id, name, i;
        node(int id, int name, int i)
            : id(id), name(name), i(i) { }
        void print() const;
    };

    struct mosfet {
        enum ElemType_t {
            nmos = 1,
            pmos = 0
        };
        class UnsupportedMOSFETType { };
        const circuit & c;
        ElemType_t ElemType;
        node *NodeD, *NodeS, *NodeG;
        std::string name;

        // V_T: threshold voltage
        // MU: mobility
        // LAMBDA: channel-width modulator
        // C_OX: oxide capacitance
        // C_J: junction capacitance
        double W, L, V_T, MU, LAMBDA, C_OX, C_J;
        mosfet(const circuit & c, ElemType_t ElemType, node *NodeD, node *NodeS, node *NodeG, std::string name,
                double W, double L, double V_T, double MU, double C_OX, double LAMBDA, double C_J)
            : c(c), ElemType(ElemType), NodeD(NodeD), NodeS(NodeS), NodeG(NodeG), name(name),
                W(W), L(L), V_T(V_T), MU(MU), C_OX(C_OX), LAMBDA(LAMBDA), C_J(C_J) { }
        void print() const;
    };

    struct linelem {
        enum ElemType_t {
            R = 'R',
            C = 'C',
            L = 'L',
            V = 'V',
            I = 'I'
        };
        class UnsupportedLinearElementType { };
        const circuit & c;
        ElemType_t ElemType;
        std::string name;
        node *Node1, *Node2;
        linelem(const circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2)
            : c(c), ElemType(ElemType), name(name), Node1(Node1), Node2(Node2) { }
        virtual void print() const;
    };
    struct resistor : linelem {
        double resistance;
        resistor(const circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2, double resistance)
            : linelem(c, ElemType, name, Node1, Node2), resistance(resistance) { }
        void print() const;
    };
    struct storage_device : linelem {
        storage_device(const circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2)
            : linelem(c, ElemType, name, Node1, Node2) { }
    };
    struct capacitor : storage_device {
        double capacitance, initial_voltage;
        capacitor(const circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2, double capacitance, double initial_voltage)
            : storage_device(c, ElemType, name, Node1, Node2), capacitance(capacitance), initial_voltage(initial_voltage) { }
        void print() const;
    };
    struct inductor : storage_device {
        double inductance, initial_current;
        inductor(const circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2, double inductance, double initial_current)
            : storage_device(c, ElemType, name, Node1, Node2), inductance(inductance), initial_current(initial_current) { }
        void print() const;
    };
    struct power_source : linelem {
        enum TYPE_t {
            DC = 0,
            // AC = 1,
            PWL = 2
        };
        class UnsupportedPowerSourceType { };
        TYPE_t SOURCE_TYPE;
        power_source(const circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2, TYPE_t SOURCE_TYPE)
            : linelem(c, ElemType, name, Node1, Node2), SOURCE_TYPE(SOURCE_TYPE) { }
        void print() const;
    };
    struct V_source : power_source {
        V_source(const circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2, TYPE_t SOURCE_TYPE)
            : power_source(c, ElemType, name, Node1, Node2, SOURCE_TYPE) { }
        void print() const;
    };
    struct V_dc : V_source {
        double voltage_value;
        V_dc(const circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2, TYPE_t SOURCE_TYPE, double voltage_value)
            : V_source(c, ElemType, name, Node1, Node2, SOURCE_TYPE), voltage_value(voltage_value) { }
    };
    struct V_pwl : V_source {
        typedef std::vector< std::pair<double, double> > voltages_t;
        voltages_t voltages;
        V_pwl(const circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2, TYPE_t SOURCE_TYPE, voltages_t voltages)
            : V_source(c, ElemType, name, Node1, Node2, SOURCE_TYPE), voltages(voltages) { }
    };
    struct I_source : power_source {
        I_source(const circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2, TYPE_t SOURCE_TYPE)
            : power_source(c, ElemType, name, Node1, Node2, SOURCE_TYPE) { }
        void print() const;
    };
    struct I_dc : I_source {
        double current_value;
        I_dc(const circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2, TYPE_t SOURCE_TYPE, double current_value)
            : I_source(c, ElemType, name, Node1, Node2, SOURCE_TYPE), current_value(current_value) { }
    };
    struct I_pwl : I_source {
        typedef std::vector< std::pair<double, double> > currents_t;
        currents_t currents;
        I_pwl(const circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2, TYPE_t SOURCE_TYPE, currents_t currents)
            : I_source(c, ElemType, name, Node1, Node2, SOURCE_TYPE), currents(currents) { }
    };


    // static variables
    static node * const gnd;


    // functions
    circuit();
    circuit(const std::string & filename);
    ~circuit();

    analysis * run(matlab * const m) const;
    void print() const;
    size_t step_num;

private:

    // elements and nodes
    std::vector<linelem*> linelems;
    std::vector<mosfet*> mosfets;
    std::unordered_map<int,node*> nodes;

    // analysis
    void * analysis_type; // analysis::TYPE_t *
    std::vector<int> PLOTNV, PLOTBV, PLOTBI;
    double time_step;
    double stop_time;
    mutable analysis * run_pointer;

    // friends
    friend class linelem;
    friend class analysis;
    friend class dc;
    friend class tran;
    friend class helper;
};

#endif
