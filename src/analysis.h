
#ifndef __ANALYSIS_H
#define __ANALYSIS_H

#include "circuit.h"
#include <vector>

class analysis {
public:
    analysis(const circuit * const c);
    enum TYPE_t {
        DC = 1, // fix
        TRAN = 3
    };
protected:
    static constexpr double precision = 0.000000001;
    const circuit * const c;
    std::vector<circuit::linelem*> itrelems;
};

class tran;
class dc;

#include "tran.h"
#include "dc.h"

#endif
