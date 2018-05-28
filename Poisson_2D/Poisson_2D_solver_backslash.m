% Testing 2D Poisson setup and solve

%Sets up the matrix equation for solving 2D Poisson equation with optional space
%varying dielectric constants. Boundary conditions are:

%-a fixed voltage at (x,0) and (x, Nz) defined by V_bottomBC and V_topBC
%which are defining the  electrodes

%-insulating boundary conditions: V(0,z) = V(1,z) and V(0,N+1) = V(1,N) (N
%is the last INTERIOR mesh point.
%so the potential at the boundary is assumed to be the same as just inside
%the boundary. Gradient of potential normal to these boundaries is 0.

%Matrix equation AV*V = bV where AV is a sparse matrix
%(generated using spdiag), V is the solution for electric potential, and bV
%is the rhs which contains the charge densities and boundary conditions

clear all

global num_cell N num_elements

%% Physical Constants
q =  1.60217646*10^-19;         %elementary charge, C
kb = 1.3806503*10^-23;          %Boltzmann const., J/k
T = 296.;                      %temperature:  Koster says they use room temperature...: I find that btw. 293 and 300 get no effect on JV curve...
epsilon_0 =  8.85418782*10^-12; %F/m
Vt = (kb*T)/q;

%% System Setup

num_cell = 10;
N = num_cell -1;   %number of INTERIOR mesh points
num_elements = N^2;  %NOTE: this will specify number of elements in the solution vector V which = (num_cell +1 -2)^2 b/c we are in 2D
%(num_cell +1) = # of mesh pts, b/c matlab starts indexing from 1, then -2
%the endpts

Va = 0.; %applied voltage

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

%% Define boundary conditions and initial conditions
V_bottomBC = 0;
V_topBC = Va/Vt;

% Initial conditions
diff = (V_topBC - V_bottomBC)/num_cell;  
index = 0;
for j = 1:N  %corresponds to z coord
    index = index +1;
    V(index) = diff*i;
    for i = 2:N  %elements along the x direction assumed to have same V
        index = index +1;
        V(index) = V(index-1);
    end
end
        
%Side BCs will be filled in from V, since are insulating BC's
%i's in these BC's correspond to the x-value (z values are along a line,
%top and bottom)
for i = 1:N
    V_leftBC(i) = V(i);  %just corresponds to values of V in the 1st subblock
    V_rightBC(i) = V(num_elements - N + i);
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