function pid = pid_broja(P)
% calculate PID with Ibroja implementation from dit
s = size(P);
Nvar = length(s);
if Nvar ~= 3
    error('only 3 variables supported')
end

thism = mfilename('fullpath');
[folder, name, ext] = fileparts(thism);

fname = [tempname '.mat'];
save(fname,'P')
cmd = sprintf('python %s/pidbroja.py %s',folder,fname);
status = system(cmd);

if status~=0
    error('python returned an error code: %d',status)
end

dat = load(fname);
pid = dat.pid;

delete(fname)
