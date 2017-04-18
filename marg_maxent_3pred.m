function [Pme] = marg_maxent_3pred(P) 
s = size(P);
Nvar = length(s);
% if Nvar ~= 3
%     error('only 3 variables supported')
% end

save pyP P
!python mme3pred.py
dat = load('pyP.mat');
Pme = dat.Pme;
