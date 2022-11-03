
#ifndef __CIRCUIT_H
#define __CIRCUIT_H

#include <string>
#include <vector>
#include <unordered_map>
#include <Eigen/Eigen>
#include <utility>

class circuit {
public:
    // structs
    struct node {
        int id, name, i;
        std::vector<double> voltages;
        node(int id, int name, int i)
            : id(id), name(name), i(i) { }
        double voltage(const int & t = -1) const;
        void print() const;
    };

    struct diode {
        circuit & c;
        node *Node1, *Node2;
        // std::string name;
        std::vector<double> voltages, currents;
        diode(circuit & c, node *Node1, node *Node2)
            : c(c), Node1(Node1), Node2(Node2) { }

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
        circuit & c;
        ElemType_t ElemType;
        std::string name;
        node *Node1, *Node2;
        linelem(circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2)
            : c(c), ElemType(ElemType), name(name), Node1(Node1), Node2(Node2) { }
        virtual double voltage(const int & t = -1) const;
        virtual void print() const;
    };
    struct resistor : linelem {
        double resistance;
        resistor(circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2, double resistance)
            : linelem(c, ElemType, name, Node1, Node2), resistance(resistance) { }
        double current(const int & t = -1) const;
        void print() const;
    };
    struct storage_device : linelem {
        double current(const int & t = -1) const;
        std::vector<double> currents;
        storage_device(circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2)
            : linelem(c, ElemType, name, Node1, Node2) { }
    };
    struct capacitor : storage_device {
        double capacitance, initial_voltage;
        capacitor(circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2, double capacitance, double initial_voltage)
            : storage_device(c, ElemType, name, Node1, Node2), capacitance(capacitance), initial_voltage(initial_voltage) { }
        double conductance() const;
        void print() const;
    };
    struct inductor : storage_device {
        double inductance, initial_current;
        inductor(circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2, double inductance, double initial_current)
            : storage_device(c, ElemType, name, Node1, Node2), inductance(inductance), initial_current(initial_current) { }
        double conductance() const;
        void print() const;
    };
    struct power_source : linelem {
        enum TYPE_t {
            DC = 0,
            // AC = 1,
            PWL = 2
        };
        TYPE_t SOURCE_TYPE;
        power_source(circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2, TYPE_t SOURCE_TYPE)
            : linelem(c, ElemType, name, Node1, Node2), SOURCE_TYPE(SOURCE_TYPE) { }
        void print() const;
    };
    struct V_source : power_source {
        virtual double voltage(const int & t = -1) const = 0;
        double current(const int & t = -1) const;
        std::vector<double> currents;
        V_source(circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2, TYPE_t SOURCE_TYPE)
            : power_source(c, ElemType, name, Node1, Node2, SOURCE_TYPE) { }
        void print() const;
    };
    struct V_dc : V_source {
        double voltage(const int & t = -1) const;
        double voltage_value;
        V_dc(circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2, TYPE_t SOURCE_TYPE, double voltage_value)
            : V_source(c, ElemType, name, Node1, Node2, SOURCE_TYPE), voltage_value(voltage_value) { }
    };
    struct V_pwl : V_source {
        double voltage(const int & t = -1) const;
        typedef std::vector< std::pair<double, double> > voltages_t;
        voltages_t voltages;
        V_pwl(circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2, TYPE_t SOURCE_TYPE, voltages_t voltages)
            : V_source(c, ElemType, name, Node1, Node2, SOURCE_TYPE), voltages(voltages) { }
    };
    struct I_source : power_source {
        virtual double current(const int & t = -1) const = 0;
        I_source(circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2, TYPE_t SOURCE_TYPE)
            : power_source(c, ElemType, name, Node1, Node2, SOURCE_TYPE) { }
        void print() const;
    };
    struct I_dc : I_source {
        double current(const int & t = -1) const;
        double current_value;
        I_dc(circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2, TYPE_t SOURCE_TYPE, double current_value)
            : I_source(c, ElemType, name, Node1, Node2, SOURCE_TYPE), current_value(current_value) { }
    };
    struct I_pwl : I_source {
        double current(const int & t = -1) const;
        typedef std::vector< std::pair<double, double> > currents_t;
        currents_t currents;
        I_pwl(circuit & c, ElemType_t ElemType, std::string name, node *Node1, node *Node2, TYPE_t SOURCE_TYPE, currents_t currents)
            : I_source(c, ElemType, name, Node1, Node2, SOURCE_TYPE), currents(currents) { }
    };
    // struct V_pwl : V_source { };
    // struct I_pwl : I_source { };


    // static variables
    static node * const gnd;


    // functions
    circuit();
    circuit(const std::string & filename);
    ~circuit();

    void dc();
    void tran();

    void print() const;
    void to_json(const std::string & filename) const;
    size_t step_num;

private:

    typedef Eigen::SparseLU< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver_t;

    // elements and nodes
    std::vector<linelem*> linelems;
    std::vector<diode*> diodes;
    std::unordered_map<int,node*> nodes;

    // tran
    std::vector<int> PLOTNV, PLOTBV, PLOTBI;
    double time_step;
    double stop_time;
    void tran_fill_A(solver_t & A, const size_t & n, const size_t & m) const;
    void tran_step(const solver_t & A, const size_t & n, const size_t & m);

    // friends
    friend class circuit_interface;
    friend class linelem;
};

#endif
