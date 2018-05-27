% Setup of 2D Poisson matrix using loop and spdiag

function AV = SetAV_2D(epsilon)

global num_elements N

AV_val = zeros(num_elements, 5);   %this is a matrix which will just store the non-zero diagonals of 2D Poisson matrix

%These are listed in order from lowest diagonal to highest diagonal

for i = 1:N*(N-1)      %this is the lowest diagonal (1st element corresponds to Nth row  (number of elements = N*(N-1)
    AV_val(i,1) = -1.;          
end

%NOTE: this is tricky!-->some elements are 0 (at the corners of the
%subblocks)
for i = 1:num_elements-1      %this is the lower diagonal (below main diagonal) (1st element corresponds to 2nd row)
    if(mod(i, N) == 0)
        AV_val(i,2) = 0;   %these are the elements at subblock corners
    else
        AV_val(i,2) = -1.;  
    end
end

for i =  1:num_elements      %main diagonal
    AV_val(i,3) = 4.;
end

for i = 2:num_elements      %main uppper diagonal, matlab fills this from the bottom (so i = 2 corresponds to 1st row in matrix)
    if(mod(i-1,N) ==0)       %i-1 b/c indexed from 2
        AV_val(i,4) = 0;
    else
        AV_val(i,4) = -1.;
    end
end

for i = 1+N:num_elements      %far upper diagonal, matlab fills from bottom, so this starts at 1+N (since 1st element is in the 2nd subblock of matrix)
    AV_val(i,5) = -1.;            %1st element corresponds to 1st row.   this has N^2 -N elements
end

%all not specified elements will remain zero, as they were initialized
%above.

AV = spdiags(AV_val, [-N -1 0 1 N], num_elements, num_elements); %A = spdiags(B,d,m,n) creates an m-by-n sparse matrix by taking the columns of B and placing them along the diagonals specified by d.
 %diagonals  [-N -1 0 1 N].  -N and N b/c the far diagonals are in next
 %subblocks, N diagonals away from main diag.
                                                                                                                          
                                                                  
                                                                   