function Pb = pid_broja(P)
% calculate PID with Ibroja implementation from dit
s = size(P);
Nvar = length(s);
if Nvar ~= 3
    error('only 3 variables supported')
end

fname = tempname;
save(fname,'P')
cmd = sprintf('python pidbrojadist.py %s',fname);
status = system(cmd);

if status~=0
    error('python returned an error code: %d',status)
end

dat = load(fname);
Pb = dat.Pb;

delete([fname '.mat'])
