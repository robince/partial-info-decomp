function P2 = marg_maxent2(P) 
s = size(P);
Nvar = length(s);

fname = tempname;
save(fname,'P')
cmd = sprintf('python mme2.py %s',fname);
status = system(cmd);

if status~=0
    error('python returned an error code: %d',status)
end

dat = load(fname);
P2 = dat.P2;

delete([fname '.mat'])