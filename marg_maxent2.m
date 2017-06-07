function P2 = marg_maxent2(P) 
s = size(P);
Nvar = length(s);

fname = tempname;
save(fname,'P')

thism = mfilename('fullpath');
[folder, name, ext] = fileparts(thism);

cmd = sprintf('python %s/mme2.py %s',folder,fname);
status = system(cmd);

if status~=0
    error('python returned an error code: %d',status)
end

dat = load(fname);
P2 = dat.P2;

delete([fname '.mat'])