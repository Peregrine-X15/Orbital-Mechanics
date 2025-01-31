function [U_o, U_1, U_2, U_3] = UnivFns_iii(alpha, chi, tol)

U_o = 1;
U_1 = chi;

dUo = U_o;
io = 1;
while abs(dUo) > tol
    dUo = ((-1)^io*(alpha*chi^2)^io)/factorial(2*io);
    U_o = U_o + dUo;
    io = io + 1;
end

dU1 = U_1;
i1 = 1;
while abs(dU1) > tol
    dU1 = chi*((-1)^i1)*(alpha*chi^2)^i1/factorial(1+2*i1);
    U_1 = U_1 + dU1;
    i1 = i1 + 1;
end


U_2 = (1-U_o)/alpha;
U_3 = (chi-U_1)/alpha;