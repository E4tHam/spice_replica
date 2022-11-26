
#include "circuit.h"
#include "matlab.h"
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

int main(int argc, char const *argv[]) {
    matlab m;
    string js = m.ckt_to_json(argv[1]);

    circuit a{js};
    // a.dc();
    a.tran();
    // a.print();
    a.to_json("out.json");
    return 0;
}
