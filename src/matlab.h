
#ifndef __MATLAB_H
#define __MATLAB_H

#include "engine.h"
#include "matrix.h"
#include <string>
#include <iostream>


class matlab {
    public:
        matlab() {
            std::cout << "Starting MATLAB engine...";
            if (!(ep = engOpen(""))) {
                std::cerr << "FAILED!" << std::endl;
                exit(1);
            } else {
                std::cout << "SUCCESS!" << std::endl;
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
