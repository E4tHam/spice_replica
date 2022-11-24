
# ESspice

<https://github.com/sifferman/ESspice>

UCSB ECE 594BB F22 by Peng Li

## About

This is a simplified open-source SPICE implementation.

### Supported Devices

* DC Voltage Sources
* PWL Voltage Sources
* DC Current Sources
* PWL Current Sources
* Resistors
* Capacitors
* Inductors
* N-MOSFETs
* P-MOSFETs

### Supported Analysis Types

* DC
* Transient (via Forward Euler Approximation)

## Running The Code Using The Provided `main.cpp`

1. In MATLAB, run the `"ckt_to_json.m"` script to parse the `".ckt"` files into `".json"` files.
2. Build the C++ code using `make`.
3. Run the C++ code with `./main <json file>` (ex. `./main circuits_json/nand3.json`) to create an `"out.json"` file.
4. In MATLAB, run the `"json_plot.m"` script to plot the specified nodes in `"out.json"`.

## To Do

* Analysis class
  * DC analysis subclass
  * Tran analysis subclass
  * DC and TRAN have different matrix solvers
* For each type of analysis, each device has
  * a "`stamp()`" function that stamps its infuence onto a solver matrix
  * a "`get_voltage()`" and "`get_current()`" function that retrieves its voltage and current from the solved matrix
* AC Analysis
