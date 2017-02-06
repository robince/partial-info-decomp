function compare_ped(P)
% utility function to compare different information decompositions

lat = lattice2d();

% Imin
lat = calc_pi_wb(lat, P, @Imin);
PImin = lat.PI;
arr = sprintf('%6.4f ',PImin);
fprintf(1,'PI with Imin    : [ %s ]\n',arr);

% Ibroja
pid = pid_broja(P);
arr = sprintf('%6.4f ',pid);
fprintf(1,'PI with Ibroja  : [ %s ]\n',arr);

% Iccs
lat = calc_pi(lat, P, @Iccs);
PIccs = lat.PI;
arr = sprintf('%6.4f ',PIccs);
fprintf(1,'PI with Iccs    : [ %s ]\n',arr);

% PED
lat = lattice3d();
lat = calc_pe(lat,P,@Hcs);

pid = pid_from_ped(lat);
arr = sprintf('%6.4f ',pid);
fprintf(1,'PI from PED     : [ %s ]\n',arr);

pid = monopid_from_ped(lat);
arr = sprintf('%6.4f ',pid);
fprintf(1,'Mono-PI from PED: [ %s ]\n',arr);

% isclosefp = @(a,b) abs(a - b) <= 1e-15;
% assert(isclosefp(sum(lat.PI),lat.Icap(end)));

fprintf(1,'\n');