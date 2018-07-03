# Quadrotor control scripts
## Test environment
* OS: Windows Server 2012 R2
* Compiler: Microsoft Visual C++ 2015 Professional
* CPU: Intel Xeon E5-2650 V4 @2.2G Hz (2 CPUs, 24 cores in total)
* Hyper-Threading and Turbo Boost disabled


## Using ACADO
1. Copy the folder `ACADO\Quadrotor` to `ACADO\interfaces\matlab\examples\code_generation\nmpc\`.
2. Run `Quadrotor.m` to define the NMPC problem.
3. Run `Quadrotor_sim.m` to do the closed-loop simulation.

## Using AutoGenU
1. Run `AutoGenU_Quadcopter.mw`.

## Using ParNMPC
1. Find the quadrotor example under the  `ParNMPC-Beta2.0` branch.
2. Follow the user manual.
