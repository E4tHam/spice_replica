
#include "circuit.h"
#include "matlab.h"
#include "analysis.h"
#include "helper.h"
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

int main(int argc, char const *argv[]) {

    if (argc!=2) {
        cerr << "Usage: ./main <.ckt file>" << endl;
        return 1;
    }

    // create matlab instance
    matlab m;

    // load circuit
    string js = m.ckt_to_json(argv[1]);
    circuit a{js};

    // dc analysis
    dc dc_a(&a);

    // tran analysis
    tran tran_a(&a, 1.0e-11, 2.0e-8);
    helper::export_circuit(&tran_a, "out.json"); // to replace

    cout << "Done." << endl;
    return 0;
}
