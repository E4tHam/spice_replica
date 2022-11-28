
#include "matlab.h"

#include <string>

std::string matlab::ckt_to_json(const std::string & filename) {

    mxArray * mxFilename = mxCreateString(&filename[0]);
    engPutVariable(ep, "filename", mxFilename);
    engEvalString(ep, "json_string = ckt_to_json(filename);");
    std::string json_string = mxArrayToString(engGetVariable(ep, "json_string"));

    return json_string;
}

void matlab::show_plot(const std::vector<double> & data, const std::string & title, const std::string & xlabel, const std::string & ylabel, const double & xtick, const double & xlim) {
    mxArray * mxData = mxCreateUninitNumericMatrix(data.size(), 1, mxDOUBLE_CLASS, mxREAL);
    for (size_t i = 0; i < data.size(); i++) {
        mxGetDoubles(mxData)[i] = data[i];
    }
    mxArray * mxTitle = mxCreateString(&title[0]);
    mxArray * mxXLabel = mxCreateString(&xlabel[0]);
    mxArray * mxYLabel = mxCreateString(&ylabel[0]);
    mxArray * mxXTick = mxCreateDoubleScalar(xtick);
    mxArray * mxXLim = mxCreateDoubleScalar(xlim);
    engPutVariable(ep, "in_data", mxData);
    engPutVariable(ep, "in_title", mxTitle);
    engPutVariable(ep, "in_xlabel", mxXLabel);
    engPutVariable(ep, "in_ylabel", mxYLabel);
    engPutVariable(ep, "in_xtick", mxXTick);
    engPutVariable(ep, "in_xlim", mxXLim);
    engEvalString(ep, "show_plot(in_data, in_title, in_xlabel, in_ylabel, in_xtick, in_xlim);");
}
