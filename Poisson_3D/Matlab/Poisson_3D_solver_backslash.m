% Testing 3D Poisson setup and solve

%Sets up the matrix equation for solving 3D Poisson equation with optional space
%varying dielectric constants. Boundary conditions are:

%-a fixed voltage at (x,y, 0) and (x, y, Nz) defined by V_bottomBC and V_topBC
%which are defining the  electrodes

%-insulating boundary conditions: V(0, y, z) = V(1, y, z) and V(x, 0, z) =  V(x, 1, z) 
%                                 V(N+1, y, z) = V(N, y, z) and V(x, N+1,
%                                 z) = V(x, N, z)

%N is the last INTERIOR mesh point.
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
tic

num_cell = 4;
N = num_cell -1;   %number of INTERIOR mesh points. (note: total # of mesh pts is N+1)
num_elements = N^3;  %NOTE: this will specify number of elements in the solution vector V which = (num_cell +1 -2)^2 b/c we are in 2D
%(num_cell +1) = # of mesh pts, b/c matlab starts indexing from 1, then -2
%the endpts

Va = 1.; %applied voltage

%Matrix of the system's net charge
%for now make the net charge = 0  --> this will give a linear 3D potential
%plot
netcharge = zeros(num_elements,1);  %make 1D matrix, since easier to incorporate into bV

%% Define dielectric constant matrix
%NOTE: epsilons will be defined at 1/2 integer points, so epsilons inside
%the cells, not at cell boundaries
%will use indexing: i + 1/2 is defined as i+1 for the index
epsilon = ones(num_cell+1, num_cell+1, num_cell+1); 
%later will fill with real epsilons corresponding to the different layers

%NOTE: I unfortunately can't define epsilon(0,..) in matlab, so the endpts
%are at 1...


%% Define boundary conditions and initial conditions
V_bottomBC = 0;  %z = 0
V_topBC = Va/Vt; % z = N+1

% Initial conditions
diff = (V_topBC - V_bottomBC)/num_cell;
index = 0;
for k = 1:N
    index = index +1;
    V(index) = diff*k;
    for cnt = 2:N^2  %elements along the x and y direction
        index = index +1;
        
        V(index) = V(index-1);
    end
end


%Side BCs will be filled in from V, since are insulating BC's
%i's in these BC's correspond to the x-value (z values are along a line,
%top and bottom)
%THESE NEED TO BE UPDATED AT EVERY ITERATION OF POISSON SOLVE
%side BC's are now matrices (2D)
index = 0;
for k = 1:N
    for j = 1:N  %1st vertical subblock set--> the j's iterate
        V_leftBC_x(j, k) = V(index + (j-1)*N + 1);  %just corresponds to values of V in the 1st subblock
        V_rightBC_x(j, k) = V(index + j*N);
    end
    index = index+N*N;  %brings us to next vertical subblock set
end

%y bc's taken from beginning and end subblocks of each vertical subblock
%set
index = 0;
for k = 1:N
    for i =  1:N  %
        V_leftBC_y(i, k) = V(index + i);  
        V_rightBC_y(i, k) = V(index + i + N*N - N);
    end
    index = index + N*N;
end


%% Set up matrix equation and solve
AV = SetAV_3D(epsilon);
%spy(AV);  %allows to see matrix structure, very useful!

%set up rhs of Poisson equation. Note for epsilons, are assuming that
%epsilons at the boundaries are the same as espilon 1 cell into interior of
%device
CV = 1;  %for now

bV = CV*netcharge;

index = 0;
for k = 1:N
    if (k == 1)  %--------------------------------------------------------------------
        for j = 1:N
            if(j == 1)  %different for 1st subblock
                for i = 1:N
                    index = index +1;
                    if (i == 1)  %1st element has 2 BC's
                        bV(index,1) = bV(index,1) + epsilon(0+1, 1+1, 1+1)*V_leftBC_x(1,1) + epsilon(1+1, 0+1, 1+1)*V_leftBC_y(1,1)  + epsilon(1+1, 1+1, 0+1)*V_bottomBC;  %+1 b/c netcharge and epsilon include endpoints but i,j index only the interior
                    elseif (i==N)
                        bV(index,1) = bV(index,1) + epsilon(N+1+1, 1+1, 1+1)*V_rightBC_x(1,1) + epsilon(N+1, 0+1, 1+1)*V_leftBC_y(N,1) + epsilon(N+1, 1+1, 0+1)*V_bottomBC;
                    else  %middle elements in 1st subblock
                        bV(index,1) = bV(index,1) + epsilon(i+1, 0+1, 1+1)*V_leftBC_y(i, 1) + epsilon(i+1, 1+1, 0+1)*V_bottomBC;
                    end
                end
            elseif(j == N)  %different for last subblock within k=1 subblock group
                for i = 1:N
                    index = index +1;
                    if (i==1)  %1st element has 2 BC's
                        bV(index,1) = bV(index,1) + epsilon(0+1,N+1, 1+1)*V_leftBC_x(N,1) + epsilon(1+1, N+1+1, 1+1)*V_rightBC_y(1,1) + epsilon(1+1, N+1, 0+1)*V_bottomBC;
                    elseif (i==N)
                        bV(index,1) = bV(index,1) + epsilon(N+1+1,N+1, 1+1)*V_rightBC_x(N,1) + epsilon(N+1, N+1+1, 1+1)*V_rightBC_y(N,1) + epsilon(N+1,N+1, 0+1)*V_bottomBC;
                    else %inner rows of Nth subblock
                        bV(index,1) = bV(index,1) + epsilon(i+1,N+1+1, 1+1)*V_rightBC_y(i,1) + epsilon(i+1, N+1, 0+1)*V_bottomBC;
                    end
                end
            else  %interior subblocks of 1st (k = 1) subblock group
                for i = 1:N
                    index = index +1;
                    if (i==1)
                        bV(index,1) = bV(index,1) + epsilon(0+1, j+1, 1+1)*V_leftBC_x(j,1) + epsilon(1+1, j+1, 0+1)*V_bottomBC;
                    elseif (i==N)
                        bV(index,1) = bV(index,1) + epsilon(N+1+1, j+1, 1+1)*V_rightBC_x(j,1) + epsilon(N+1, j+1, 0+1)*V_bottomBC;
                    else
                        bV(index,1) = bV(index,1) + epsilon(i+1, j+1, 0+1)*V_bottomBC;
                    end
                end
            end
        end
    elseif (k == N)  %last subblock group
        for j = 1:N
            if(j == 1)  %different for 1st subblock
                for i = 1:N
                    index = index +1;
                    if (i == 1)  %1st element has 2 BC's
                        bV(index,1) = bV(index,1) + epsilon(0+1, 1+1, N+1)*V_leftBC_x(1,N) + epsilon(1+1, 0+1, N+1)*V_leftBC_y(1,N)  + epsilon(1+1, 1+1, N+1+1)*V_topBC;  %+1 b/c netcharge and epsilon include endpoints but i,j index only the interior
                    elseif (i==N)
                        bV(index,1) = bV(index,1) + epsilon(N+1+1, 1+1, N+1)*V_rightBC_x(1,N) + epsilon(N+1, 0+1, N+1)*V_leftBC_y(N,N) + epsilon(N+1, 1+1, N+1+1)*V_topBC;
                    else  %middle elements in 1st subblock
                        bV(index,1) = bV(index,1) + epsilon(i+1, 0+1, N+1)*V_leftBC_y(i, N) + epsilon(i+1, 1+1, N+1+1)*V_topBC;
                    end
                end
            elseif(j == N)  %different for last subblock within k=N subblock group
                for i = 1:N
                    index = index +1;
                    if (i==1)  %1st element has 2 BC's
                        bV(index,1) = bV(index,1) + epsilon(0+1,N+1, N+1)*V_leftBC_x(N,N) + epsilon(1+1, N+1+1, N+1)*V_rightBC_y(1,N) + epsilon(1+1, N+1, N+1+1)*V_topBC;
                    elseif (i==N)
                        bV(index,1) = bV(index,1) + epsilon(N+1+1,N+1, N+1)*V_rightBC_x(N,N) + epsilon(N+1, N+1+1, N+1)*V_rightBC_y(N,N) + epsilon(N+1,N+1, N+1+1)*V_topBC;
                    else %inner rows of Nth subblock
                        bV(index,1) = bV(index,1) + epsilon(i+1,N+1+1, N+1)*V_rightBC_y(i,N) + epsilon(i+1, N+1, N+1+1)*V_topBC;
                    end
                end
            else  %interior subblocks of last (k = N) subblock group
                for i = 1:N
                    index = index +1;
                    if (i==1)
                        bV(index,1) = bV(index,1) + epsilon(0+1, j+1, N+1)*V_leftBC_x(j,N) + epsilon(1+1, j+1, N+1+1)*V_topBC;
                    elseif (i==N)
                        bV(index,1) = bV(index,1) + epsilon(N+1+1, j+1, N+1)*V_rightBC_x(j,N) + epsilon(N+1, j+1, N+1+1)*V_topBC;
                    else
                        bV(index,1) = bV(index,1) + epsilon(i+1, j+1, N+1+1)*V_topBC;
                    end
                end
            end
        end
    else   %interior subblock groups (k=2:N-1)
        for j = 1:N
            if(j == 1)  %different for 1st subblock
                for i = 1:N
                    index = index +1;
                    if (i == 1)  %1st element has 2 BC's
                        bV(index,1) = bV(index,1) + epsilon(0+1, 1+1, k+1)*V_leftBC_x(1,k) + epsilon(1+1, 0+1, k+1)*V_leftBC_y(1,k);  %+1 b/c netcharge and epsilon include endpoints but i,j index only the interior
                    elseif (i==N)
                        bV(index,1) = bV(index,1) + epsilon(N+1+1, 1+1, k+1)*V_rightBC_x(1,k) + epsilon(N+1, 0+1, k+1)*V_leftBC_y(N,k);
                    else  %middle elements in 1st subblock
                        bV(index,1) = bV(index,1) + epsilon(i+1, 0+1, k+1)*V_leftBC_y(i, k);
                    end
                end
            elseif(j == N)  %different for last subblock
                for i = 1:N
                    index = index +1;
                    if (i==1)  %1st element has 2 BC's
                        bV(index,1) = bV(index,1) + epsilon(0+1,N+1, k+1)*V_leftBC_x(N,k) + epsilon(1+1, N+1+1, k+1)*V_rightBC_y(1,k);
                    elseif (i==N)
                        bV(index,1) = bV(index,1) + epsilon(N+1+1,N+1, k+1)*V_rightBC_x(N,k) + epsilon(N+1, N+1+1, k+1)*V_rightBC_y(N,k);
                    else %inner rows of Nth subblock
                        bV(index,1) = bV(index,1) + epsilon(i+1,N+1+1, k+1)*V_rightBC_y(i,k);
                    end
                end
            else  %interior subblocks
                for i = 1:N
                    index = index +1;
                    if (i==1)
                        bV(index,1) = bV(index,1) + epsilon(0+1, j+1, k+1)*V_leftBC_x(j,k);
                    elseif (i==N)
                        bV(index,1) = bV(index,1) + epsilon(N+1+1, j+1, k+1)*V_rightBC_x(j,k);
                    end
                end
            end
        end
    end
end
        
                
          

%solve for V
%spparms('spumoni',2)
V = AV\bV;

toc

%plot V
% surf(1:N,1:N,reshape(V,N,N))  %need to reshape the vector into a 2D matrix, to make a 2D plot
