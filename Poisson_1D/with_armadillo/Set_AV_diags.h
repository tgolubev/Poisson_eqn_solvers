#ifndef SET_AV_DIAGS_H
#define SET_AV_DIAGS_H

#include<vector>


std::vector<double> set_main_Vdiag(std::vector<double> &epsilon, std::vector<double> &a);
std::vector<double> set_upper_Vdiag(std::vector<double> &epsilon, std::vector<double> &b);
std::vector<double> set_lower_Vdiag(std::vector<double> &epsilon, std::vector<double> &c);

#endif // SET_AV_DIAGS_H
