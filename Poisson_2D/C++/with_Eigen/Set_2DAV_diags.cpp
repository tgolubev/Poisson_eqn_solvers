#include <iostream>
#include<vector>

#include "Set_2DAV_diags.h"

//far lower diagonal
std::vector<double> set_far_lower_Vdiag(Eigen::MatrixXd &epsilon, std::vector<double> &Farlower_diag){
    for(int index = 1; index<= N*(N-1); index++){
        int i = index % N;
        if(i==0) i=N;
        int j = 1+floor((index-1)/N);

        Farlower_diag[index] = -(epsilon(i,j) + epsilon(i+1,j))/2.;
    }
    return Farlower_diag;
}



//lower diag
std::vector<double> set_lower_Vdiag(Eigen::MatrixXd &epsilon, std::vector<double> &lower_diag){
    //int num_elements = lower_diag.size() -1;

    for (int index = 1; index<=num_elements-1;index++){  //      %this is the lower diagonal (below main diagonal) (1st element corresponds to 2nd row)
        int i = index % N;        // %this is x index of V which element corresponds to (note if this = 0, means these are the elements which are 0);
        int j = 1 + floor((index-1)/N);

        if(index % N == 0)
            lower_diag[index] = 0; //  %these are the elements at subblock corners
        else
            lower_diag[index] = -(epsilon(i,j) + epsilon(i,j+1))/2.;
    }

    return lower_diag;

}


//main diagonal
std::vector<double> set_main_Vdiag(Eigen::MatrixXd &epsilon, std::vector<double> &main_diag){
    //int num_elements = main_diag.size() - 1;

    for (int index =  1; index<=num_elements;index++){//      %main diagonal
        int i = index % N;
        if(i ==0)        //        %the multiples of N correspond to last index
            i = N;
        int j = 1 + floor((index-1)/N);

        main_diag[index] = (epsilon(i+1,j) + epsilon(i+1,j+1))/2. + (epsilon(i,j) + epsilon(i,j+1))/2. + (epsilon(i,j+1) + epsilon(i+1,j+1))/2. + (epsilon(i,j) + epsilon(i+1,j))/2.;
    }

    return main_diag;
}

//upper diagonal
std::vector<double> set_upper_Vdiag(Eigen::MatrixXd &epsilon, std::vector<double> &upper_diag){
    //int num_elements = upper_diag.size() -1;  //-1 since vectors fill from 0, but I'm indexing from 1
    //NOTE: num_elements here, is actually # of elements in the off-diagonal
    
    for (int index = 1; index <=num_elements-1;index++) {  //      %main uppper diagonal, matlab fills this from the bottom (so i = 2 corresponds to 1st row in matrix)
        int i = index % N;
        int j = 1 + floor((index-1)/N);

        if(index  % N ==0)
            upper_diag[index] = 0;
        else
            upper_diag[index] = -(epsilon(i+1,j) + epsilon(i+1,j+1))/2.;
   }

    return upper_diag;
}


//far upper diagonal
std::vector<double> set_far_upper_Vdiag(Eigen::MatrixXd &epsilon, std::vector<double> &far_upper_diag){

    for (int index = 1; index <= num_elements-N; index++) { //      %far upper diagonal, matlab fills from bottom, so this starts at 1+N (since 1st element is in the 2nd subblock of matrix)
         int i = index % N;
        if(i ==0)      //          %the multiples of N correspond to last index
            i = N;
        int j = 1 + floor((index-N)/N);

         far_upper_diag[index] = -(epsilon(i,j+1) + epsilon(i+1,j+1))/2.;        //    %1st element corresponds to 1st row.   this has N^2 -N elements
    }
    return far_upper_diag;

}

