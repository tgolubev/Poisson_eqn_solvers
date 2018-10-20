#include <iostream>
#include<vector>


//this is a in tridiag_solver
std::vector<double> set_main_Vdiag(std::vector<double> &epsilon, std::vector<double> &main_diag){
    int num_elements = main_diag.size() - 1;

    for (int i=1; i<= num_elements;i++){
        main_diag[i] = -2.*epsilon[i];
    }
    return main_diag;
}

//this is b in tridiag_solver
std::vector<double> set_upper_Vdiag(std::vector<double> &epsilon, std::vector<double> &upper_diag){
    int num_elements = upper_diag.size() -1;  //-1 since vectors fill from 0, but I'm indexing from 1
    //NOTE: num_elements here, is actually # of elements in the off-diagonal
    
    for (int i = 1; i<=num_elements; i++){
        upper_diag[i] = epsilon[i];
    }
    return upper_diag;
}

//this is c in tridiag_solver
std::vector<double> set_lower_Vdiag(std::vector<double> &epsilon, std::vector<double> &lower_diag){
    int num_elements = lower_diag.size() -1;
    
    for (int i = 1; i<=num_elements; i++){
        lower_diag[i] = epsilon[i];
    }
    return lower_diag;
    
}
