function lat = calc_pi_ri(lat,Pjoint,Icap)
% Calculate PI on a redundancy lattice using alternative sumation method
% that takes into account the problem of additivity for non-disjoint child
% nodes

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

% PI at each node
lat.PI = NaN(size(lat.Icap));

% ascend through levels of the lattice
Nlevels = max(lat.level);
for li=1:Nlevels
    nodes = find(lat.level==li);
    % reset backprop factors for children of this level
    lat.backprop = ones(size(lat.Icap));
    for ni=nodes
        lat = calc_pi_node(lat,ni);
    end
    % apply backprop factor for children of this level
    lat.PI = lat.PI.*lat.backprop;
end


function lat = calc_pi_node(lat,ni)
% incoming information deltas for this node
deltas = zeros(size(lat.children{ni}));
if isempty(deltas)
    % no children
    lat.PI(ni) = lat.Icap(ni);
    return
end
for ci=1:length(deltas)
    deltas(ci) = lat.Icap(ni) - lat.PI(lat.children{ni}(ci));
end
% ignore negative deltas (??)
deltas(deltas<0) = 0;
pi = min(deltas);
all_children = unique(recurse_children(lat,ni,[]));
if pi>0
    % contributing input nodes
    idx = isclosefp(deltas,pi);
    % check if inputs are disjoint
    % in higher order lattices might need to do this more carefully
    % with 3d lattice every node with multiple inputs has all inputs
    % non-disjoint
    % with 2d lattice every node with multiple inputs (only {12}) has all
    % inputs disjoint
    % number of non-disjoint inputs
    N = sum(idx);
    input_nodes = lat.children{ni}(idx);
    if lat.Nx==3
        Nnondis = N;
        nondis_input_nodes = input_nodes;
    elseif lat.Nx==2
        Nnondis = 1;
        nondis_input_nodes = [];
    end
    

    if N>1
        % normalise pi
        % should it be by N inputs or sum inputs??
        % how should combine factors from multiple parents?
        pi = (N*pi) ./ Nnondis;
        lat.backprop(nondis_input_nodes) = lat.backprop(nondis_input_nodes) .* (1./Nnondis);
    end
    all_children = setdiff(all_children,nondis_input_nodes);
end

thsPI = pi - sum(lat.PI(all_children));
lat.PI(ni) = max(thsPI,0);
% ni
% lat.PI
% if ni==4
%     keyboard
% end

function children = recurse_children(lat,ni,children)
children = [children lat.children{ni}];
for ci=lat.children{ni}
    children = recurse_children(lat,ci,children);
end

function y = isclosefp(a,b)
y = abs(a - b) <= eps(max(abs(a), abs(b)));