
#include "circuit.h"

int main(int argc, char const *argv[]) {
    circuit a{argv[1]};
    a.print();
    return 0;
}
