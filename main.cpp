
#include "circuit.h"
#include <iostream>
#include <iomanip>
#include <string>

int main(int argc, char const *argv[]) {
    // std::cout << std::fixed << std::setprecision(2);
    circuit a{argv[1]};
    a.print();
    for (int i = 0; i < std::stoi(argv[2]); i++) {
        a.step();
        a.print();
    }
    a.to_json("out.json");
    return 0;
}
