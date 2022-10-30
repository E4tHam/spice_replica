
#ifndef __CIRCUIT_INTERFACE_H
#define __CIRCUIT_INTERFACE_H

#include "circuit.h"

struct circuit_interface {
    static void circuit_from_filename(circuit * c, const std::string & filename);
    static void export_circuit(const circuit * const c, const std::string & filename);
};

#endif
