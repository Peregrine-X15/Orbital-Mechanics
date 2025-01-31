function [E, Ehis, errhis] = solvekepler(Mt, e)

%global e P;
%M = 2*pi*t/P;

E0 = Mt;
err = E0 - e*sin(E0) - Mt;
thres = 1e-9;
Ecurr = Mt;
iter = 1;
MAXiter = 200;
Ehis = E0;
errhis = err;
correction = 1;
while abs(correction) > thres
    correction = (Mt - (Ecurr - e*sin(Ecurr)))/(1 - e*cos(Ecurr));
    Enew = Ecurr + correction;
    err = Mt - (Enew - e*sin(Enew));
    Ecurr = Enew;
    Ehis = [Ehis; Enew];
    errhis = [errhis; err];
    iter = iter + 1;
    if iter > MAXiter
        break;
    end
end

E = Ecurr;