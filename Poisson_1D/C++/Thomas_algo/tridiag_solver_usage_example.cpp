//This is a code for testing the tridiagonal solver in the context of solving the 1D Poisson equation
//Allows to test both TriCRSSolver and Thomas solver


//March 30, 2018
//Timofey Golubev

#include <iostream>
#include "Set_AV_diags.h"
#include <vector>
#include "thomas_tridiag_solve.h"
#include <chrono>

// valid indices in a vector fall into the range 0 to vectorSize - 1. So if want to fill from 1 to num_elements, need to make a num_elements + 1 sized vector


int main()
{
    const int num_elements = 9999;  //number of nodes in the system
    std::vector<double> epsilon(num_elements+2);
    std::vector<double> a (num_elements+1);//main diag
    std::vector<double> b(num_elements); //upper diag, size can be = num_elements b/c off-diags are 1 element less than  main diag
    std::vector<double> c(num_elements);//lower diag
    std::vector<double> V(num_elements+1); //vector for solution, electric potential
    std::vector<double> rhs(num_elements+1);

    double V_leftBC = 0.;
    double V_rightBC= 0.9;

    for(int i=0;i<=num_elements+1;i++){
        epsilon[i] = 3.8;
    }

    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();  //start clock timer

    a = set_main_Vdiag(epsilon, a);
    b = set_upper_Vdiag(epsilon, b);
    c = set_lower_Vdiag(epsilon, c);

    //setup rhs of Poisson eqn.
    for (int i = 1;i<= num_elements;i++){
         rhs[i] = 0;  //for now
    }

    //test i.e. having a dipole at accross 24-25th node
    rhs[24] = 0.681; //corresponds to having 5 holes/(75nm)^2 plane at the interface, with z-mesh size of 1nm
    rhs[25] = -0.681;

    //for bndrys
    rhs[1] = rhs[1] - epsilon[1]*V_leftBC;
    rhs[num_elements] = rhs[num_elements] - epsilon[num_elements]*V_rightBC;

    //V = TriCRSSolver(a, b, c, rhs);
     V = Thomas_solve(a, b, c, rhs);

    std::chrono::high_resolution_clock::time_point finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> time = std::chrono::duration_cast<std::chrono::duration<double>>(finish-start);

    //output the result to terminal

    for(int i = 1;i<= num_elements; i++){
        std::cout << V[i] << " " <<std::endl;
    }

    std::cout << "Total CPU time = " << time.count() << std::endl;


    return 0;
}
