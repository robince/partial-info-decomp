function lat = calc_pi_mvn(lat,Cfull,varsizes,Icap,forcenn,normlevels)
% Calculate PI on a redundancy lattice using Williams and Beer summation
% inputs:
% lat - lattice structure
% Cfull - covariance matrix of MVN system
% varsizes - split of dimensions of full covariance into predictor and
% target (target always ordered last). ie [2 2 1] means 2 2d predictors, 1d
% target (size(Cfull)=[5 5])
% forcenz - threshold negative values on the lattice (default false)
% normlevels - normalize values across levels of the lattice (default
% false)
%
% if only lat provided calculate PI using existing Icap
% otherwise, recalculate Icap

if nargin>1
    if sum(varsizes) ~= size(Cfull,1)
        error('wrong number of variables specified')
    end
    if length(varsizes)~=3 && length(varsizes)~=4
        error('only 2 or 3 predictors supported')
    end

    % calc Icap for each node
    for ni=1:lat.Nnodes
        lat.Icap(ni) = Icap(lat.A{ni}, Cfull, varsizes);
    end
end

if lat.Nx>3
    error('calc_pi: too many variables')
end

if nargin<5
    forcenn = false;
end
if nargin<6
    normlevels = false;
end

% use equation (7) from Williams and Beer to calculate
% PI at each node
lat.PI = NaN(size(lat.Icap));
% raw PI before non-disjoint normalisation
lat.PIraw = NaN(size(lat.Icap));

% ascend through levels of the lattice
Nlevels = max(lat.level);
for li=1:(Nlevels-1)
    nodes = find(lat.level==li);
    for ni=nodes
        lat = calc_pi_node(lat,ni,forcenn,normlevels);
    end
end
% don't enforce non-negativitity for top node
lat = calc_pi_node(lat,lat.top,false,normlevels);


function lat = calc_pi_node(lat,ni,nonneg,normlevels)
if nargin<3
    nonneg = false;
end
children = lat.children{ni};
if isempty(children)
    % no children
    thsPI = lat.Icap(ni);
    if nonneg
        thsPI = max(thsPI,0);
    end
    lat.PI(ni) = thsPI;
    lat.PIraw(ni) = thsPI;
    return
end
all_children = recurse_children(lat,ni,[]);

if normlevels
    PIchildren = normalise_levels(lat, all_children);
else
    PIchildren = lat.PI(all_children);
end
thsPI = lat.Icap(ni) - sum(PIchildren);
if nonneg
    thsPI = max(thsPI,0);
end

lat.PI(ni) = thsPI;
lat.PIraw(ni) = thsPI;

if ni==lat.top
    lat.PI(all_children) = PIchildren;
end



function normPI = normalise_levels(lat,children)
% normalise to correct for non-additivity of non-disjoint nodes

% values for this set of children
PIraw = lat.PIraw(children);
levels = lat.level(children);
labels = lat.labels(children);
A = lat.A(children);
normPI = PIraw;

for li=1:lat.Nlevels
    nodes = find(levels==li);
    levelPI = PIraw(nodes);
    posPInodes = nodes(abs(levelPI)>1e-15);
    posPIvars = A(posPInodes);
    posPIvars = cell2mat([posPIvars{:}]);
    if length(posPIvars) ~= length(unique(posPIvars))
        % have non-disjoint positive PI contributions at this level
        
        % using structure of 3rd order lattice (might need more logic to
        % determine pairwise disjoint-ness for higher order lattices)
        if li==4
            % special case level 4 for 3 variable lattice
            % one node contains all variables
            fullnode = find(strcmpi(labels,'{12}{13}{23}'));
            if isempty(fullnode) || PIraw(fullnode)<1e-15
                % all sources at this level are disjoint so no
                % normalisation required
                continue
            elseif length(posPInodes)==1
                % only {12}{13}{23} is non-zero so no normalization
                % required
                continue
            end
            % only normalise by 2 here even if more posPInodes, because
            % there are only 2 disjoint copies at this level
            normPI(posPInodes) = PIraw(posPInodes) ./ 2;
        else
            normPI(posPInodes) = PIraw(posPInodes) ./ length(posPInodes);
        end
    end
end


function children = recurse_children(lat,ni,children)
children = [children lat.children{ni}];
for ci=lat.children{ni}
    children = recurse_children(lat,ci,children);
end
children = unique(children);

