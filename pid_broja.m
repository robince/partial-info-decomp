function pid = pid_broja(P)
% calculate PID with Ibroja implementation from dit
s = size(P);
Nvar = length(s);
if Nvar ~= 3
    error('only 3 variables supported')
end

save pyP P
!python pidbroja.py
dat = load('pyP.mat');
pid = dat.pid;