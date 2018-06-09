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
T = 296.;                      %temperature
epsilon_0 =  8.85418782*10^-12; %F/m
Vt = (kb*T)/q;

%% System Setup

num_cell = 10;
N = num_cell -1;   %number of INTERIOR mesh points
num_elements = N^2;  %NOTE: this will specify number of elements in the solution vector V which = (num_cell +1 -2)^2 b/c we are in 2D
%(num_cell +1) = # of mesh pts, b/c matlab starts indexing from 1, then -2
%the endpts

Va = 1.; %applied voltage

%Matrix of the system's net charge
%for now make the net charge = 0  --> this will give a linear 2D potential
%plot
netcharge = zeros(num_elements,num_elements);
%netcharge = ones(num_elements,num_elements);  %this will make a slightly curved 2D potential plot

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
    V(index) = diff*j;
    for i = 2:N  %elements along the x direction assumed to have same V
        index = index +1;
        
        V(index) = V(index-1);
    end
end

%Side BCs will be filled in from V, since are insulating BC's
%i's in these BC's correspond to the x-value (z values are along a line,
%top and bottom)
%THESE NEED TO BE UPDATED AT EVERY ITERATION OF POISSON SOLVE
for i = 1:N
    V_leftBC(i) = V((i-1)*N + 1);  %just corresponds to values of V in the 1st subblock
    V_rightBC(i) = V(i*N);
end


%% Set up matrix equation and solve
AV = SetAV_2D(epsilon);
%spy(AV);  %allows to see matrix structure, very useful!

%set up rhs of Poisson equation. Note for epsilons, are assuming that
%epsilons at the boundaries are the same as espilon 1 cell into interior of
%device
index = 0;
for j = 1:N
    if(j ==1)  %different for 1st subblock
        for i = 1:N
            index = index +1;
            if (i==1)  %1st element has 2 BC's
                bV(index,1) = netcharge(i,j) + epsilon(i,j)*(V_leftBC(1) + V_bottomBC);
            elseif (i==N)
                bV(index,1) = netcharge(i,j) + epsilon(i,j)*(V_rightBC(1) + V_bottomBC);
            else
                bV(index,1) = netcharge(i,j) + epsilon(i,j)*V_bottomBC;
            end
        end
    elseif(j == N)  %different for last subblock
        for i = 1:N
            index = index +1;
            if (i==1)  %1st element has 2 BC's
                bV(index,1) = netcharge(i,j) + epsilon(i,j)*(V_leftBC(N) + V_topBC);
            elseif (i==N)
                bV(index,1) = netcharge(i,j) + epsilon(i,j)*(V_rightBC(N) + V_topBC);
            else
                bV(index,1) = netcharge(i,j) + epsilon(i,j)*V_topBC;
            end
        end
    else %interior subblocks
        for i = 1:N
            index = index +1;
            if(i==1)
                bV(index,1) = netcharge(i,j) + epsilon(i,j)*V_leftBC(j);
            elseif(i==N)
                bV(index,1) = netcharge(i,j) + epsilon(i,j)*V_rightBC(j);
            else
                bV(index,1) = netcharge(i,j);
            end
        end
    end
end

%solve for V
V = AV\bV

%plot V
surf(1:N,1:N,reshape(V,N,N))  %need to reshape the vector into a 2D matrix, to make a 2D plot
