//This is a code for solving the 2D Poisson equation with possibly position dependent dielectric
//constant using Armadillo (linear algebra library)

//June 2018
//Timofey Golubev

#include <iostream>
#include "Set_2DAV_diags.h"
#include <vector>

#include "parameters.h"

#include "armadillo.h"  //note: armadillo features are located in namespace arma


int main()
{
    double Va = 1.;

    //arma vectors and matrices fill from 0.
    arma::sp_mat AV(num_elements, num_elements);   //create an armadillo sparse matrix
    arma::vec vec_main_diag (num_elements);
    arma::vec vec_upper_diag(num_elements-1);
    arma::vec vec_lower_diag(num_elements-1);
    arma::vec vec_far_lower_diag(num_elements-N);
    arma::vec vec_far_upper_diag(num_elements-N);
    arma::vec vec_rhs(num_elements);
    arma::vec vec_V(num_elements);

    arma::mat epsilon(num_cell+2,num_cell+2);
    epsilon.fill(1.); //set all elements = 1;
    arma::mat netcharge(num_cell+2,num_cell+2);
    netcharge.fill(0.);  //for now make net charge = 0

    std::vector<double> main_diag (num_elements+1);
    std::vector<double> upper_diag(num_elements);
    std::vector<double> lower_diag(num_elements);
    std::vector<double> far_lower_diag(num_elements-N+1);
    std::vector<double> far_upper_diag(num_elements-N+1);

    std::vector<double> V(num_elements+1); //vector for solution, electric potential
    std::vector<double> rhs(num_elements+1);
    //these are still filled as before (using from i = 1)

    //Define boundary conditions and initial conditions
    double V_bottomBC = 0.;
    double V_topBC= Va/Vt;

    //Initial conditions
    double diff = (V_topBC - V_bottomBC)/num_cell;
    int index = 0;
    for (int j = 1;j<=N;j++){//  %corresponds to z coord
        index++;
        V[index] = diff*j;
        for (int i = 2;i<=N;i++){//  %elements along the x direction assumed to have same V
            index++;
            V[index] = V[index-1];
        }
    }

    //side BCs, insulating BC's
    std::vector<double> V_leftBC(N+1);
    std::vector<double> V_rightBC(N+1);

    for(int i = 1;i<=N;i++){
        V_leftBC[i] = V[(i-1)*N +  1];
        V_rightBC[i] = V[i*N];
    }

    //Set up matrix equation
    main_diag = set_main_Vdiag(epsilon, main_diag);
    upper_diag = set_upper_Vdiag(epsilon, upper_diag);
    lower_diag = set_lower_Vdiag(epsilon, lower_diag);
    far_upper_diag = set_far_upper_Vdiag(epsilon, far_upper_diag);
    far_lower_diag = set_far_lower_Vdiag(epsilon, far_lower_diag);

    //setup rhs of Poisson eqn.
    int index2 = 0;
    for(int j = 1;j<=N;j++){
        if(j==1){
            for(int i = 1;i<=N;i++){
                index2++;                //THIS COULD BE MODIFIES TO a switch ,case statement--> might be cleaner
                if(i==1){
                    rhs[index2] = netcharge(i,j) + epsilon(i,j)*(V_leftBC[1] + V_bottomBC);
                }else if(i == N)
                    rhs[index2] = netcharge(i,j) + epsilon(i,j)*(V_rightBC[1] + V_bottomBC);
                else
                    rhs[index2] = netcharge(i,j) + epsilon(i,j)*V_bottomBC;
            }
        }else if(j==N){
            for(int i = 1; i<=N;i++){
                index2++;
                if(i==1)
                    rhs[index2] = netcharge(i,j) + epsilon(i,j)*(V_leftBC[N] + V_topBC);
                else if(i == N)
                    rhs[index2] = netcharge(i,j) + epsilon(i,j)*(V_rightBC[N] + V_topBC);
                else
                    rhs[index2] = netcharge(i,j) + epsilon(i,j)*V_topBC;
            }
        }else{
            for(int i = 1;i<=N;i++){
                index2++;
                if(i==1)
                    rhs[index2] = netcharge(i,j) + epsilon(i,j)*V_leftBC[j];
                else if(i == N)
                    rhs[index2] = netcharge(i,j) + epsilon(i,j)*V_rightBC[j];
                else
                    rhs[index2] = netcharge(i,j);
            }
        }
    }

    //armadillo solving setup
    //Note: I fill my a,b,c std::vectors from 1 (more consistent with math eqns)
    //and arma fills from 0.
    for(int i = 1; i< main_diag.size(); i++){
        vec_main_diag(i-1) = main_diag[i];
        vec_rhs(i-1) = rhs[i];
    }
    std::cout << vec_rhs <<std::endl;
    for(int i = 1;i< upper_diag.size();i++){
        vec_upper_diag(i-1) = upper_diag[i];
        vec_lower_diag(i-1) = lower_diag[i];

    }
    for(int i = 1;i< far_upper_diag.size();i++){
        vec_far_upper_diag(i-1) = far_upper_diag[i];
        vec_far_lower_diag(i-1) = far_lower_diag[i];
    }


    //fill the diagonals of AV using .diag(k). Argument specifies the diagonal (k>0 are upper diags, and k<0 are lower diags).
    //Note: rhs must  be an arma vec for diag to work.
    AV.diag(0) = vec_main_diag;
    AV.diag(1) = vec_upper_diag;
    AV.diag(-1) = vec_lower_diag;
    AV.diag(N) = vec_far_upper_diag;
    AV.diag(-N) = vec_far_lower_diag;

    std::cout << AV << std::endl;  //armadillo has overloaded << to output matrices
    //std::cout<< vec_a << std::endl;
    //std::cout<< vec_b << std::endl;
    //std::cout<< vec_c << std::endl;
   //std::cout<< vec_rhs << std::endl;

    vec_V = spsolve(AV, vec_rhs, "lapack");

    //output the result to terminal
    std:: cout << "solution " << vec_V << std::endl;


    return 0;
}
