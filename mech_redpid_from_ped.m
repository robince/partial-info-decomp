function pid = pid_from_ped(lat)
% calculate PID of redundant entropy information measure

if lat.Nx ~= 3
    error('only 2 variable PID supported')
end

pid = zeros(1,5);

n = lat.nodes;
mech = abs( min( lat.PI(n('{1}{2}')), 0) );
% source redundancy
pid(1) = lat.PI(n('{1}{2}{3}')) - mech;
% mechanistic redundancy
pid(2) = mech;

% unique info
pid(3) = lat.PI(n('{1}{3}'));
   
pid(4) = lat.PI(n('{2}{3}'));

% synergy
pid(5) = lat.PI(n('{3}{12}'));