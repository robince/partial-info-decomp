function lat = lattice2d()
% build 2d redundancy lattice

lat.Nx = 2;
lat.Nnodes = 4;
lat.Nelements = 3;
lat.elements = {1, 2, [1 2]};

lat.A = { {1 2}, {1}, {2}, {[1 2]} };
% lat.labels = { '{1}{2}', '{1}', '{2}', '{12}' }
lat.level = [1 2 2 3];
labels = cell(size(lat.A));
for ni=1:lat.Nnodes
    thsA = lat.A{ni};
    st = cell(size(thsA));
    for ai=1:length(st)
        st{ai} = sprintf('{%s}',sprintf('%d',thsA{ai}));
    end
    labels{ni} = sprintf('%s',st{:});
end
lat.labels = labels;

nodes = containers.Map(lat.labels, 1:lat.Nnodes);
parents = cell(1,lat.Nnodes);
parents{nodes('{1}{2}')} = [nodes('{1}') nodes('{2}')];
parents{nodes('{1}')} = nodes('{12}');
parents{nodes('{2}')} = nodes('{12}');

% build children for birectional link
children = cell(1,lat.Nnodes);
for ni=1:lat.Nnodes
    for parent=parents{ni}
        children{parent} = [children{parent} ni];
    end
end

lat.nodes = nodes;
lat.parents = parents;
lat.children = children;

lat.top = nodes('{12}');
lat.bottom = nodes('{1}{2}');
lat.Icap = NaN(1,lat.Nnodes);
lat.PI = NaN(1,lat.Nnodes);

