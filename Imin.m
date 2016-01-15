function Imin = Imin(A, Pjoint)

% A - cell array of sources
% Pjoint - full joint distribution

s = size(Pjoint);
Sm = s(end); % number of target values
Nx = length(s)-1; % number of dependent variables
vars = 1:Nx;

NA = length(A);
Pele(NA).Pa = []; % intialize struct

% Ps
Ps = Pjoint;
for xi=1:Nx
    Ps = squeeze(sum(Ps,1));
end

% build distributions for each source
for ai=1:NA
    thsA = A{ai};
    % vars to sum over
    sumover = setdiff(vars, thsA);
    Pas = Pjoint;
    for ii=1:length(sumover)
        Pas = sum(Pas, sumover(ii));
    end
    Pas = squeeze(Pas);
    thsPs = reshape(Ps,[ones(1, length(size(Pas))-1) Sm]);
    Pele(ai).PaCs = bsxfun(@rdivide, Pas, thsPs);
    Pele(ai).Pa = squeeze(sum(Pas,length(size(Pas))));
end

% calculate Imin
Dkl = @(x,y) sum(x(y~=0).*(log2(x(y~=0))-log2(y(y~=0))));

% loop over target values
minIs = zeros(1,Sm);
for si=1:Sm
    Is = zeros(1,NA);
    for ai=1:NA
        % deal with matlabs horrible indexing contraints
        % make s axis first
        s = length(size(Pele(ai).PaCs));
        PaCs = permute(Pele(ai).PaCs, [s 1:s-1]);
        Is(ai) = Dkl(PaCs(si,:)', Pele(ai).Pa(:));
    end
    minIs(si) = min(Is);
end
Imin = sum(Ps.*minIs);