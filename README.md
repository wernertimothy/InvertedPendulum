# Inverted Pendulum
simulation of different controlers for the inverted pendulum stabilization.

## Dependencies
To run the code, [qpOASES](https://www.coin-or.org/qpOASES/doc/3.2/doxygen/index.html) isrequired.
Download it [here](https://github.com/coin-or/qpOASES), unzip and inside the directory .../interfaces/matlab run the command in Matlab

´´´matlab
make
´´´

This will create two mex-files which will allow you to call the qpOASES solver.
To see further information see the [Manual](https://www.coin-or.org/qpOASES/doc/3.0/manual.pdf).

## Simulation
Currently running simulaitons are an LQR and a linear tracking MPC.

## ToDo's
* add noise and estimation
* do the swing-up