for i = 1:num_elements
    diagonal(i) = -2*epsilon(i+1)/epsilon_0;
end
 diagonal(l_HTL_int) = -(epsilon(l_HTL_int+1)/epsilon_0 + epsilon(l_HTL_int)/epsilon_0);
 diagonal(l_ETL_int) = -(epsilon(l_ETL_int+1)/epsilon_0 + epsilon(l_ETL_int)/epsilon_0);
 diagonal(l_BCP_int)= -(epsilon(l_BCP_int+1)/epsilon_0 + epsilon(l_BCP_int)/epsilon_0);
    
    
   for i = 1:num_elements 
     b(i) = epsilon(i+1)/epsilon_0;
   end
     b(l_HTL_int) = epsilon(l_HTL_int -1)/epsilon_0;% //takes the value to the left  is . . * case
    b(l_ETL_int)= epsilon(l_ETL_int-1)/epsilon_0;
    b(l_BCP_int) = epsilon(l_BCP_int-1)/epsilon_0;
    for i = 1:num_elements  
        c(i) =epsilon(i+1)/epsilon_0;
    end
    
rhs = bV;
%//Forward substition
for  i = 2:num_elements
    cdiag_ratio = c(i-1)/diagonal(i-1);  %i-1 b/c need to use c1 and a1, and i starts at 2...
    diagonal(i) = diagonal(i) - cdiag_ratio*b(i-1);
    rhs(i) = rhs(i) - cdiag_ratio*rhs(i-1);
end

%//Backward substitution
Vsolution(num_elements) = rhs(num_elements)/diagonal(num_elements); %//linear eqn corresponding to last row
for i = num_elements:-1:2 Vsolution(i-1) = (rhs(i-1)-Vsolution(i)*b(i-1))/diagonal(i-1);
end

Vsolution