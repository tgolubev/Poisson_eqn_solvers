//This is a code for testing tridiagonal matrix solving using Armadillo (linear algebra library)

//June 2018
//Timofey Golubev

#include <iostream>
#include "Set_AV_diags.h"
#include <vector>

#include "armadillo.h"  //note: armadillo features are located in namespace arma

// valid indices in a vector fall into the range 0 to vectorSize - 1. So if want to fill from 1 to num_cell, need to make a num_cell + 1 sized vector

int main()
{

    const int num_cell = 100;  //number of nodes in the system

    //arma vectors and matrices fill from 0.
    arma::sp_mat AV(num_cell-1, num_cell-1);   //create an armadillo sparse matrix
    arma::vec vec_a(num_cell-1);
    arma::vec vec_b(num_cell-2);
    arma::vec vec_c(num_cell-2);
    arma::vec vec_rhs(num_cell-1);
    arma::vec vec_V(num_cell-1);

    std::vector<double> epsilon(num_cell+1);
    std::vector<double> a (num_cell);//main diag
    std::vector<double> b(num_cell-1); //upper diag, size can be = num_cell b/c off-diags are 1 element less than  main diag
    std::vector<double> c(num_cell-1);//lower diag
    std::vector<double> V(num_cell); //vector for solution, electric potential
    std::vector<double> rhs(num_cell);
    //these are still filled as before (using from i = 1)

    double V_leftBC = 0.;
    double V_rightBC= 0.9;

    for(int i=0;i<=num_cell;i++){
        epsilon[i] = 3.8;
    }

    a = set_main_Vdiag(epsilon, a);
    b = set_upper_Vdiag(epsilon, b);
    c = set_lower_Vdiag(epsilon, c);

    //setup rhs of Poisson eqn.
    for (int i = 1;i<= num_cell;i++){
         rhs[i] = 0;  //for now
    }
    //test i.e. having a dipole at accross 24-25th node
    rhs[24] = 0.681; //corresponds to having 5 holes/(75nm)^2 plane at the interface, with z-mesh size of 1nm
    rhs[25] = -0.681;

    //for bndrys
    rhs[1] = rhs[1] - epsilon[0]*V_leftBC;
    rhs[num_cell-1] = rhs[num_cell-1] - epsilon[num_cell]*V_rightBC;

    //armadillo solving setup
    //Note: I fill my a,b,c std::vectors from 1 (more consistent with math eqns)
    //and arma fills from 0.
    for(int i = 1; i< a.size(); i++){
        vec_a(i-1) = a[i];
        vec_rhs(i-1) = rhs[i];
    }
    for(int i = 1;i< b.size();i++){
        vec_b(i-1) = b[i];
        vec_c(i-1) = c[i];

    }

    //fill the diagonals of AV using .diag(k). Argument specifies the diagonal (k>0 are upper diags, and k<0 are lower diags).
    //Note: rhs must  be an arma vec for diag to work.
    AV.diag(0) = vec_a;
    AV.diag(1) = vec_b;
    AV.diag(-1) = vec_c;

    std::cout << AV << std::endl;  //armadillo has overloaded << to output matrices
    std::cout<< vec_a << std::endl;
    std::cout<< vec_b << std::endl;
    std::cout<< vec_c << std::endl;
    std::cout<< vec_rhs << std::endl;

    vec_V = spsolve(AV, vec_rhs, "lapack");

    //V = TriCRSSolver(a, b, c, rhs);
    //V = Thomas_solve(a, b, c, rhs);

    //output the result to terminal
    std:: cout << "solution " << vec_V << std::endl;

    //this loop for cout, only works for std::vector. It crashes for arma vectors
    /*
    for(int i = 0;i<= num_cell-2; i++){
        std::cout << V[i] << " " <<std::endl;
    }
    */


    return 0;
}
