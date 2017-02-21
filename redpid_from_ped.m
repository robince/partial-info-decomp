function pid = pid_from_ped(lat)
% calculate PID of redundant entropy information measure

if lat.Nx ~= 3
    error('only 2 variable PID supported')
end

pid = zeros(1,4);
n = lat.nodes;

% redundancy
pid(1) = lat.PI(n('{1}{2}{3}'));
   
% unique info
pid(2) = lat.PI(n('{1}{3}'));
   
pid(3) = lat.PI(n('{2}{3}'));

% synergy
pid(4) = lat.PI(n('{3}{12}'));