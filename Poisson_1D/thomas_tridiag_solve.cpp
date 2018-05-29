#include <vector>

// diagonal = array containing elements of main diagonal. indices: (a1.....an)
// b = array containing elements of upper diagonal. indices (b1....b_n-1)
//c = array containing elements of lower diagonal. indices (c1...c_n-1)
std::vector<double> Thomas_solve(std::vector<double> &diagonal, std::vector<double> &b, std::vector<double> &c, std::vector<double> &rhs){

//This uses Thomas algorithm for tridiagonal matrix (special case of Gaussian elimination)
int num_elements = diagonal.size()-1;  //matrices indexed from 0, but we use it from 1 here
double cdiag_ratio;
std::vector<double> x(num_elements+2);

//Forward substition
for (int i = 2; i <= num_elements; i++) {
    cdiag_ratio = c[i-1]/diagonal[i-1];   //i-1 b/c need to use c1 and a1, and i starts at 2...
    diagonal[i] -= cdiag_ratio*b[i-1];
    rhs[i] -= cdiag_ratio*rhs[i-1];
}

//Backward substitution
x[num_elements] = rhs[num_elements]/diagonal[num_elements]; //linear eqn corresponding to last row
for (int i = num_elements; i >1; i--) x[i-1] = (rhs[i-1]-x[i]*b[i-1])/diagonal[i-1];

return x;
}
