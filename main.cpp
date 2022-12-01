
#include "circuit.h"
#include "matlab.h"
#include "analysis.h"
#include <iostream>
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
    cout << "Finished loading circuit." << endl;

    // Run circuit
    analysis * run = a.run(&m);

    // Allow time to look at figures
    cout << "Press return to exit...";
    getchar();
    cout << "Exiting." << endl;
    return 0;
}
