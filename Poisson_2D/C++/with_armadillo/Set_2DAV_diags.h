#ifndef SET_AV_DIAGS_H
#define SET_AV_DIAGS_H

#include<vector>
#include "armadillo.h"

#include "parameters.h"

std::vector<double> set_far_lower_Vdiag(arma::mat &epsilon, std::vector<double> &Farlower_diag);
std::vector<double> set_lower_Vdiag(arma::mat &epsilon, std::vector<double> &lower_diag);
std::vector<double> set_main_Vdiag(arma::mat &epsilon, std::vector<double> &main_diag);
std::vector<double> set_upper_Vdiag(arma::mat &epsilon, std::vector<double> &upper_diag);
std::vector<double> set_far_upper_Vdiag(arma::mat &epsilon, std::vector<double> &far_upper_diag);



#endif // SET_AV_DIAGS_H
