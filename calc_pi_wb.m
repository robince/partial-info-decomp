function lat = calc_PI(lat,Pjoint,Icap)
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
% start at the bottom
ni = lat.bottom;
lat.PI(ni) = lat.Icap(ni);
% recurse up the tree
lat = recurse_calc_parents(lat,ni);



function lat = recurse_calc_parents(lat,ni)
% calculate all parents of this node
lat = calc_parents(lat,ni);
% recurse up
for pi=lat.parents{ni}
    lat = recurse_calc_parents(lat,pi);
end



function lat = calc_parents(lat,ni)
for pi=lat.parents{ni}
    all_children = unique(recurse_children(lat,pi,[]));
    thsPI = lat.Icap(pi) - sum(lat.PI(all_children));
    lat.PI(pi) = thsPI;
end



function children = recurse_children(lat,ni,children)
children = [children lat.children{ni}];
for ci=lat.children{ni}
    children = recurse_children(lat,ci,children);
end

