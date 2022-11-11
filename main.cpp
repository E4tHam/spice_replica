
#include "circuit.h"
#include <iostream>
#include <iomanip>
#include <string>

int main(int argc, char const *argv[]) {
    circuit a{argv[1]};
    // a.dc();
    a.tran();
    // a.print();
    a.to_json("out.json");
    return 0;
}
