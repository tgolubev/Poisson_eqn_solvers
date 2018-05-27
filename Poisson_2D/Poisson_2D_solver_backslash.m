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

%% Define dielectric constant matrix
%NOTE: epsilons will be defined at 1/2 integer points, so epsilons inside
%the cells, not at cell boundaries
%will use indexing: i + 1/2 is defined as i+1 for the index
epsilon = zeros(num_cell+2, num_cell +2);
%later will fill with real epsilons corresponding to the different layers

%NOTE: I unfortunately can't define epsilon(0,..) in matlab, so the endpts
%are at 1...
for i = 1:num_cell+2
    epsilon(i,:) = 1;
end


%% Set up matrix equation and solve
AV = SetAV_2D(epsilon);
%spy(AV);  %allows to see matrix structure, very useful!

%set up rhs of Poisson equation
for i = 1:num_elements
    bV(i,1) = 1;
end

%solve for V
V = AV\bV