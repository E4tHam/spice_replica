
#include "engine.h"
#include "matrix.h"
#include <string>

std::string call_ckt_to_json(Engine * const ep, const std::string & filename) {

    mxArray * mxFilename = mxCreateString(&filename[0]);
    engPutVariable(ep, "filename", mxFilename);
    engEvalString(ep, "addpath('src/matlab');");
    engEvalString(ep, "json_string = ckt_to_json(filename);");
    std::string json_string = mxArrayToString(engGetVariable(ep, "json_string"));

    return json_string;
}
