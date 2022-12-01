
#include <vector>
#include <string>
#include <unordered_map>
struct analysis { typedef int TYPE_t; };
struct matlab;

class circuit {
public:
    struct node {
        int id, name, i;
    };

    struct mosfet {
        enum ElemType_t {
            nmos = 1,
            pmos = 0
        };
        circuit & c;
        ElemType_t ElemType;
        node *NodeD, *NodeS, *NodeG;
        std::string name;

        double W, L, V_T, MU, LAMBDA, C_OX, C_J;
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
    };
    struct resistor : linelem {
        double resistance;
    };
    struct storage_device : linelem {
    };
    struct capacitor : storage_device {
        double capacitance, initial_voltage;
    };
    struct inductor : storage_device {
        double inductance, initial_current;
    };
    struct power_source : linelem {
        enum TYPE_t {
            DC = 0,
            PWL = 2
        };
        TYPE_t SOURCE_TYPE;
    };
    struct V_source : power_source {
    };
    struct V_dc : V_source {
        double voltage_value;
    };
    struct V_pwl : V_source {
        typedef std::vector< std::pair<double, double> > voltages_t;
        voltages_t voltages;
    };
    struct I_source : power_source {
    };
    struct I_dc : I_source {
        double current_value;
    };
    struct I_pwl : I_source {
        typedef std::vector< std::pair<double, double> > currents_t;
        currents_t currents;
    };

    static node * gnd;

    circuit(std::string js);

    analysis * run(matlab * m) const;

private:
    // elements and nodes
    std::vector<linelem*> linelems;
    std::vector<mosfet*> mosfets;
    std::unordered_map<int,node*> nodes;

    // analysis
    analysis::TYPE_t * analysis_type;
    std::vector<int> PLOTNV, PLOTBV, PLOTBI;
    double time_step;
    double stop_time;
    analysis * run_pointer;
};
