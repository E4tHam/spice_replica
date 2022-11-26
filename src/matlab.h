
#ifndef __MATLAB_H
#define __MATLAB_H

#include "engine.h"
#include "matrix.h"
#include <string>
#include <iostream>


class matlab {
    public:
        matlab() {
            if (!(ep = engOpen(""))) {
                std::cerr << "Can't start MATLAB engine" << std::endl;
                exit(1);
            }
        }
        ~matlab() {
            engClose(ep);
        }
        std::string ckt_to_json(const std::string & filename);
    private:
        Engine * ep;
};

#endif
