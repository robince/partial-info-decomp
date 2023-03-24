function [pid P] = calcpid(SA, Na, SB, Nb, T, Nt, Ntrl)
% Sample probablity distributions from discrete data and use them 
% to calculate a PID with Iccs
% SA - source A (integer values in range [0, Na-1]
% Na - number of bins for SA
% SB - source B (integer values in range [0, Nb-1]
% Nb - number of bins for SB
% T  - target (integer values in range [0, Nt-1]
% Ntrl - number of trials

SA = SA(:);
SB = SB(:);
T = T(:);

if length(SA)~=Ntrl || length(SB)~=Ntrl || length(T)~=Ntrl
    error('calcpid: wrong number of trials in data')
end
if any(SA<0) || any(SA>Na-1)
    error('calcpid: problem with SA input')
end
if any(SB<0) || any(SB>Nb-1)
    error('calcpid: problem with SB input')
end
if any(T<0) || any(T>Nt-1)
    error('calcpid: problem with T input')
end

Sjoint = T*(Nb*Na) + SB*Na + SA;
Nj = Na*Nb*Nt;

P = prob(Sjoint, Nj);
P = reshape(P, [Na Nb Nt]);

% compute PID
lat = calc_pi(lattice2d(),P,@Iccs);
pid = lat.PI;
