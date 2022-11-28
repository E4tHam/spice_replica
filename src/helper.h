
#ifndef __HELPER_H
#define __HELPER_H

#include <nlohmann/json.hpp>
#include <vector>

class tran;

struct helper {
    static double pwl_value(const double & x, std::vector< std::pair<double,double> > v);

    static double bidirectional_access(const std::vector<double> & v, const int & bi);
};

#include "analysis.h"

#endif
