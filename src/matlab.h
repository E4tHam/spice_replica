
#ifndef __MATLAB_H
#define __MATLAB_H

#include "engine.h"
#include "matrix.h"
#include <string>
#include <iostream>
#include <vector>


class matlab {
    public:
        class CouldNotStartMatlabEngine { };
        matlab() {
            std::cerr << "Starting MATLAB engine...";
            if (!(ep = engOpen(""))) {
                std::cerr << "FAILED!" << std::endl;
                throw CouldNotStartMatlabEngine();
            } else {
                std::cerr << "SUCCESS!" << std::endl;
            }
            engEvalString(ep, "addpath('src/matlab');");
        }
        ~matlab() {
            engClose(ep);
        }
        std::string ckt_to_json(const std::string & filename);
        void show_plot(const std::vector<double> & data, const std::string & title, const std::string & xlabel, const std::string & ylabel, const double & xtick, const double & xlim);
    private:
        Engine * ep;
};

#endif
