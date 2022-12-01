
#include "engine.h"
#include "matrix.h"
#include <string>
#include <iostream>
#include <vector>


class matlab {
    public:
        matlab(); // start MATLAB engine
        ~matlab(); // stop MATLAB engine
        std::string ckt_to_json(std::string filename);
        void show_plot(std::vector<double> data, std::string title, std::string xlabel, std::string ylabel, double xtick, double xlim);
    private:
        Engine * ep;
};
