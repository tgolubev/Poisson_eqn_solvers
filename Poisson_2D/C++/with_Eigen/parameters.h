#ifndef PARAMETERS_H
#define PARAMETERS_H

const int num_cell = 301;
const int N = num_cell-1;
const int num_elements = N*N;

// Physical Constants
const double q =  1.60217646e-19;       //  %elementary charge, C
const double kb = 1.3806503e-23;      //    %Boltzmann const., J/k
const double T = 296.;                      //%temperature:
const double epsilon_0 =  8.85418782e-12; //%F/m
const double Vt = (kb*T)/q;



#endif // PARAMETERS_H
