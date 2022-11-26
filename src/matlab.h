
#ifndef __MATLAB_H
#define __MATLAB_H

#include "engine.h"
#include "matrix.h"
#include <string>

std::string call_ckt_to_json(Engine * const ep, const std::string & filename);

#endif
