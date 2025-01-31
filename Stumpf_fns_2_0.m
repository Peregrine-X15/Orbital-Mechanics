function [C, S]=Stumpf_fns_2_0(z,tol)

   %initiate
C = 1/2;
S = 1/6;

updateC = -z/factorial(4);
updateS = -z/factorial(5);

C = C + updateC;
S = S + updateS;


%compute C(z)
k = 2;
while abs(updateC) > tol
    updateC = (-1)^k*z.^k/factorial(2*(k + 1));
    C = C + updateC;
    k = k + 1;
end
%compute S(z)
k = 2;
while abs(updateS) > tol
    updateS = (-1)^k*z.^k/factorial(1 + 2*(k + 1));
    S = S + updateS;
    k = k + 1;
end