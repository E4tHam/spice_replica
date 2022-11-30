
#ifndef __ANALYSIS_H
#define __ANALYSIS_H

#include "circuit.h"
#include "matlab.h"
#include <vector>

class analysis {
public:
    analysis(const circuit * const c);
    enum TYPE_t {
        DC      = 0,
        TRAN_FE = 1,
        TRAN_BE = 2,
        TRAN_TR = 3
    };
    class UnsupportedAnalysisType { };
    virtual void plotnv(matlab * const m, const int & node_name) const = 0;
    virtual void printnv(const int & node_name) const = 0;
    // virtual void plotbv() const = 0;
    // virtual void plotbi() const = 0;
protected:
    // how accurate the simulation should be
    static constexpr double precision = 0.000000001;
    // Resistor added to every node
    static constexpr double R__n = (0.5/precision);
    const circuit * const c;
    std::vector<circuit::linelem*> itrelems;
};

class tran;
class dc;

#include "tran.h"
#include "dc.h"

#endif
