function lat = calc_Icap(lat, Pjoint, Icap)
% calculate Icap function (provided as funciton handle)
% at each node of the lattice for the joint distribution given
% in Pjoint

s = size(Pjoint);
if lat.Nx ~= (length(s)-1)
    error('Pjoint does not match lattice structure')
end

for ni=1:lat.Nnodes
    lat.Icap(ni) = Icap(lat.A{ni}, Pjoint);
end