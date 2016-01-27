function lat = lattice3d()
% build 3d redundancy lattice

lat.Nx = 3;
lat.Nnodes = 18;
lat.Nelements = 7;
lat.elements = {1, 2, 3, [1 2], [1 3], [2 3], [1 2 3]};

lat.A = {
    {1 2 3};
    {1 2};           {1 3};           {2 3};
    {1 [2 3]};       {2 [1 3]};       {3 [1 2]};
    {1};             {2};             {3};         {[1 2] [1 3] [2 3]};
    {[1 2] [1 3]};   {[1 2] [2 3]};   {[1 3] [2 3]};
    {[1 2]};         {[1 3]};         {[2 3]};
    {[1 2 3]} }';
lat.level = [1 2 2 2 3 3 3 4 4 4 4 5 5 5 6 6 6 7];
lat.Nlevels = max(lat.level);

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
n = nodes;
parents{n('{1}{2}{3}')} = [n('{1}{2}') n('{1}{3}') n('{2}{3}')];

parents{n('{1}{2}')} = [n('{1}{23}') n('{2}{13}')];
parents{n('{1}{3}')} = [n('{1}{23}') n('{3}{12}')];
parents{n('{2}{3}')} = [n('{2}{13}') n('{3}{12}')];

parents{n('{1}{23}')} = [n('{1}') n('{12}{13}{23}')];
parents{n('{2}{13}')} = [n('{2}') n('{12}{13}{23}')];
parents{n('{3}{12}')} = [n('{3}') n('{12}{13}{23}')];

parents{n('{1}')} = n('{12}{13}');
parents{n('{2}')} = n('{12}{23}');
parents{n('{3}')} = n('{13}{23}');
parents{n('{12}{13}{23}')} = [n('{12}{13}') n('{12}{23}') n('{13}{23}')];

parents{n('{12}{13}')} = [n('{12}') n('{13}')];
parents{n('{12}{23}')} = [n('{12}') n('{23}')];
parents{n('{13}{23}')} = [n('{13}') n('{23}')];

parents{n('{12}')} = n('{123}');
parents{n('{13}')} = n('{123}');
parents{n('{23}')} = n('{123}');

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

lat.top = nodes('{123}');
lat.bottom = nodes('{1}{2}{3}');
lat.Icap = NaN(1,lat.Nnodes);
lat.PI = NaN(1,lat.Nnodes);

