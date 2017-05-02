function lat = compare(lat,P)

lat = calc_pi(lat, P, @Imin);
PImin = lat.PI;
arr = sprintf('%6.4f ',PImin);
fprintf(1,'PI with Imin: [ %s ]\n',arr);
if lat.Nx==3
    display_nonzero_terms(lat)
end

% lat = calc_pi_ri(lat, P, @Iccs);
% PIccs = lat.PI;
% arr = sprintf('%6.4f ',PIccs);
% fprintf(1,'PI with Iccs: [ %s ]\n',arr);
% if lat.Nx==3
%     display_nonzero_terms(lat)
% end

% lat = calc_pi(lat, P, @Iccs);
% PIccs = lat.PI;
% arr = sprintf('%6.4f ',PIccs);
% fprintf(1,'PI with Iccs_Pind: [ %s ]\n',arr);
% if lat.Nx==3
%     display_nonzero_terms(lat)
% end

% 
% lat = calc_pi(lat, P, @Iccs_Pind_sepmarg, true, true);
% PIccs = lat.PI;
% arr = sprintf('%6.4f ',PIccs);
% fprintf(1,'PI with Iccs_Pind: [ %s ]\n',arr);
% if lat.Nx==3
%     display_nonzero_terms(lat)
% end

% 
% lat = calc_pi_wb(lat, P, @Iccs_pme);
% PIccs_pme = lat.PI;
% arr = sprintf('%6.4f ',PIccs_pme);
% fprintf(1,'PI with Iccs_pme: [ %s ]\n',arr);
% if lat.Nx==3
%     display_nonzero_terms(lat)
% end

lat = calc_pi(lat, P, @Iccs);
PIccs_pme = lat.PI;
arr = sprintf('%6.4f ',PIccs_pme);
fprintf(1,'PI with Iccs: [ %s ]\n',arr);
if lat.Nx==3
    display_nonzero_terms(lat)
end

isclosefp = @(a,b) abs(a - b) <= 1e-15;
assert(isclosefp(sum(lat.PI),lat.Icap(end)));

if lat.Nx<3
    pid = pid_broja(P);
    arr = sprintf('%6.4f ',pid);
    fprintf(1,'PI with Ibroja: [ %s ]\n',arr);
end
fprintf(1,'\n');