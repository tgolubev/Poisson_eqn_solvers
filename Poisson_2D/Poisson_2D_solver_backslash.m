% Testing 2D Poisson setup and solve

%Solves the matrix equation AV*V = bV where AV is a sparse matrix
%(generated using spdiag), V is the solution for electric potential, and bV
%is the rhs which contains the charge densities

clear all

global num_cell N num_elements

num_cell = 10;
N = num_cell+1 - 2;  %number of mesh points (excluding endpts)
num_elements = N^2;  %NOTE: this will specify number of elements in the solution vector V which = (num_cell +1 -2)^2 b/c we are in 2D
%(num_cell +1) = # of mesh pts, b/c matlab starts indexing from 1, then -2
%the endpts

epsilon(1:num_elements) = 1;

%Set up AV
AV = SetAV_2D(epsilon);
%spy(AV);  %allows to see matrix structure, very useful!

%set up rhs of Poisson equation
for i = 1:num_elements
    bV(i,1) = 1;
end

%solve for V
V = AV\bV