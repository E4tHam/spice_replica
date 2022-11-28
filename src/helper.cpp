
#include "helper.h"
#include <fstream>
#include <iostream>

using namespace std;

double helper::pwl_value(const double & x, std::vector< std::pair<double,double> > v) {
    if (x < 0)
        return v.front().second;
    if (x > v.back().first)
        return v.back().second;
    for (size_t i = 1; i < v.size(); i++) {
        if (x < v.at(i).first)
            return v.at(i-1).second;
    }
    return v.back().second;
}

double helper::bidirectional_access(const std::vector<double> & v, const int & bi) {
    size_t i = (bi<0) ? (bi + v.size()) : (bi);
    if (i >= v.size()) return NAN;
    return v.at(i);
}
