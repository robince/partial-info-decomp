function lat = calc_pi_wb(lat,Pjoint,Icap)
% Calculate PI on a redundancy lattice using Williams and Beer summation

% if only lat provided calculate PI using existing Icap
% otherwise, recalculate Icap
if nargin>1
    s = size(Pjoint);
    if lat.Nx ~= (length(s)-1)
        error('Pjoint does not match lattice structure')
    end

    % calc Icap for each node
    for ni=1:lat.Nnodes
        lat.Icap(ni) = Icap(lat.A{ni}, Pjoint);
    end
end

% use equation (7) from Williams and Beer to calculate
% PI at each node
lat.PI = NaN(size(lat.Icap));

% ascend through levels of the lattice
Nlevels = max(lat.level);
for li=1:Nlevels
    nodes = find(lat.level==li);
    for ni=nodes
        lat = calc_pi_node(lat,ni);
    end
end



function lat = calc_pi_node(lat,ni)
all_children = recurse_children(lat,ni,[]);
if isempty(all_children)
    thsPI = lat.Icap(ni);
else
    thsPI = lat.Icap(ni) - sum(lat.PI(all_children));
end
lat.PI(ni) = thsPI;



function children = recurse_children(lat,ni,children)
children = [children lat.children{ni}];
for ci=lat.children{ni}
    children = recurse_children(lat,ci,children);
end
children = unique(children);

