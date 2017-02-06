function [P1,P2] = marg_maxent2(P) 
s = size(P);
Nvar = length(s);
if Nvar ~= 3
    error('only 3 variables supported')
end

save pyP P
!python mme2.py
dat = load('pyP.mat');
P1 = dat.P1;
P2 = dat.P2;