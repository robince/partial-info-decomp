function [Pme] = marg_maxent_3pred(P) 
s = size(P);
Nvar = length(s);
% if Nvar ~= 3
%     error('only 3 variables supported')
% end

thism = mfilename('fullpath');
[folder, name, ext] = fileparts(thism);

fname = [tempname '.mat'];
save(fname,'P')
cmd = sprintf('python %s/mme3pred.py %s',folder,fname);
status = system(cmd);

if status~=0
    error('python returned an error code: %d',status)
end

dat = load(fname);
Pme = dat.Pme;

delete(fname)
