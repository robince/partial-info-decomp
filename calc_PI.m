function lat = calc_PI(lat)

if any(isnan(lat.Icap))
    error('Icap must be calcualted for each node')
end
lat.PI = NaN(size(lat.Icap));
% lat.cumPI = NaN(size(lat.Icap));

% use equation (7) from Williams and Beer
% so start at the bottom of the lattice

% level 1
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

