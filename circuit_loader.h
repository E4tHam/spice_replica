
#ifndef __CIRCUIT_LOADER_H
#define __CIRCUIT_LOADER_H

#include "circuit.h"

struct circuit_loader {
    static void circuit_from_filename(circuit * c, const std::string & filename);
};

#endif
