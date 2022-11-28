
#include "circuit.h"
#include "matlab.h"
#include "analysis.h"
#include "helper.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdio>

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

    // tran analysis
    tran tran_a(&a, 0.1, 3);
    tran_a.plotnv(&m, 1);

    // Allow time to look at figures
    cout << "Press return to exit..." << endl;
    getchar();
    cout << "Exiting." << endl;
    return 0;
}
