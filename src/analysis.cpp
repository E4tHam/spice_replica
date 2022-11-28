
#include "analysis.h"

analysis::analysis(const circuit * const c) : c(c) {
    for (auto m : c->mosfets) {
        itrelems.push_back(new circuit::capacitor(
            *c,
            circuit::linelem::C,
            (m->name+"__C_GS"),
            m->NodeG,
            m->NodeS,
            (0.5 * m->W * m->L * m->C_OX),
            NAN
        ));
        itrelems.push_back(new circuit::capacitor(
            *c,
            circuit::linelem::C,
            (m->name+"__C_GD"),
            m->NodeG,
            m->NodeD,
            (0.5 * m->W * m->L * m->C_OX),
            NAN
        ));
        itrelems.push_back(new circuit::capacitor(
            *c,
            circuit::linelem::C,
            (m->name+"__C_D0"),
            m->NodeD,
            circuit::gnd,
            m->C_J,
            NAN
        ));
        itrelems.push_back(new circuit::capacitor(
            *c,
            circuit::linelem::C,
            (m->name+"__C_S0"),
            m->NodeS,
            circuit::gnd,
            m->C_J,
            NAN
        ));
    }
}
