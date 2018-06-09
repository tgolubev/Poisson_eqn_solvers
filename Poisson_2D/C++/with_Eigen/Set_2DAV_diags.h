#ifndef SET_AV_DIAGS_H
#define SET_AV_DIAGS_H

#include<vector>
#include <Eigen/Dense>
#include<Eigen/Sparse>

#include<Eigen/IterativeLinearSolvers>
#include "parameters.h"

std::vector<double> set_far_lower_Vdiag(Eigen::MatrixXd &epsilon, std::vector<double> &Farlower_diag);
std::vector<double> set_lower_Vdiag(Eigen::MatrixXd &epsilon, std::vector<double> &lower_diag);
std::vector<double> set_main_Vdiag(Eigen::MatrixXd &epsilon, std::vector<double> &main_diag);
std::vector<double> set_upper_Vdiag(Eigen::MatrixXd &epsilon, std::vector<double> &upper_diag);
std::vector<double> set_far_upper_Vdiag(Eigen::MatrixXd &epsilon, std::vector<double> &far_upper_diag);



#endif // SET_AV_DIAGS_H
