
# ESspice

<https://github.com/sifferman/ESspice>

[![ESspice Demonstration](https://img.youtube.com/vi/4hBipFRTRAM/0.jpg)](https://www.youtube.com/watch?v=4hBipFRTRAM "ESspice Demonstration")

## About

This is a simplified open-source SPICE implementation completed for UCSB ECE 594BB F22 by Peng Li.

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

## Dependencies

* Linux OS
* MATLAB 2022b

Add the following to your `"~/.bashrc"`:

```bash
export MATLABROOT=/usr/local/MATLAB/R2022b
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MATLABROOT/extern/bin/glnxa64:$MATLABROOT/sys/os/glnxa64:$MATLABROOT/bin/glnxa64
```
