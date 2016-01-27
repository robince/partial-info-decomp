function lat = calc_pi_ri(lat,Pjoint,Icap)
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

if lat.Nx>3
    error('calc_pi_ri: too many variables')
end

% use equation (7) from Williams and Beer to calculate
% PI at each node
lat.PI = NaN(size(lat.Icap));
% raw PI before non-disjoint normalisation
lat.PIraw = NaN(size(lat.Icap));

% ascend through levels of the lattice
Nlevels = max(lat.level);
for li=1:Nlevels
    nodes = find(lat.level==li);
    for ni=nodes
        lat = calc_pi_node(lat,ni);
    end
    % normalise to correct for non-additivity of non-disjoint nodes
    levelPI = lat.PI(nodes);
    posPInodes = nodes(levelPI>0);
    posPIelems = lat.A(posPInodes);
    posPIelems = cell2mat([posPIelems{:}]);
    if length(posPIelems) ~= length(unique(posPIelems))
        % have non-disjoint positive PI contributions
        
        % using structure of 3rd order lattice (might need more logic to
        % determine pairwise disjoint-ness for higher order lattices)
        if li==4
            % special case level 4 for 3 variable lattice
            % one node contains all variables
            fullnode = lat.nodes('{12}{13}{23}');
            % other nodes at this level are disjoint
            if lat.PI(fullnode)>0
                % we have non-disjoint contribution at this level
                disjoint_nodes = setdiff(posPInodes, fullnode);
                disjointPI = lat.PI(disjoint_nodes);
                lat.PI(disjoint_nodes) = disjointPI .* sum(disjointPI) ./ sum(levelPI);
            end
        else
            % all sources at this level are non-disjoint
            lat.PI(posPInodes) = lat.PI(posPInodes).^2 ./ sum(levelPI);
        end
    end
end



function lat = calc_pi_node(lat,ni)
children = lat.children{ni};
if isempty(children)
    % no children
    lat.PI(ni) = lat.Icap(ni);
    lat.PIraw(ni) = lat.Icap(ni);
    return
end
all_children = recurse_children(lat,ni,[]);
if length(children)==1
    % use the raw PI for the single child node, as it is not double counted
    % here
    thsPI = lat.Icap(ni) - lat.PIraw(children) - sum(lat.PI(setdiff(all_children,children)));
else
    % use corrected PI values
    thsPI = lat.Icap(ni) - sum(lat.PI(all_children));
end
lat.PI(ni) = max(thsPI,0);
lat.PIraw(ni) = max(thsPI,0);



function children = recurse_children(lat,ni,children)
children = [children lat.children{ni}];
for ci=lat.children{ni}
    children = recurse_children(lat,ci,children);
end
children = unique(children);

