
# Spice Replica

<https://github.com/E4tHam/spice_replica>

UCSB ECE 594BB F22 by Peng Li

## To Do

* Analysis class
  * DC analysis subclass
  * Tran analysis subclass
  * DC and TRAN have different matrix solvers
* For each type of analysis, each device has
  * a "`stamp()`" function that stamps its infuence onto a solver matrix
  * a "`get_voltage()`" and "`get_current()`" function that retrieves its voltage and current from the solved matrix
* AC Analysis
