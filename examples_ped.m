% Examples for manuscript
% "The Partial Entropy Decomposition: Decomposing multivariate entropy via
% pointwise common surprisal"

%% Dyadic PED
Pxyz = zeros(4,4,4);
p = 1/8;
Pxyz(1,1,1) = p;
Pxyz(1,3,2) = p;
Pxyz(2,1,3) = p;
Pxyz(2,3,4) = p;
Pxyz(3,2,1) = p;
Pxyz(3,4,2) = p;
Pxyz(4,2,3) = p;
Pxyz(4,4,4) = p;


lat = lattice3d();
lat = calc_pe(lat,Pxyz,@Hcs);

fprintf(1,'James and Crutchfield Dyadic PED\n')
lat.PI

%% Triadic PED
Pxyz = zeros(4,4,4);
p = 1/8;
Pxyz(1,1,1) = p;
Pxyz(2,2,2) = p;
Pxyz(1,3,3) = p;
Pxyz(2,4,4) = p;
Pxyz(3,1,3) = p;
Pxyz(4,2,4) = p;
Pxyz(3,3,1) = p;
Pxyz(4,4,2) = p;

lat = lattice3d();
lat = calc_pe(lat,Pxyz,@Hcs);

fprintf(1,'James and Crutchfield Triadic PED\n')
lat.PI

%% XOR PED
% Griffith & Koch (2014) Fig. 6.5
Pxxy = zeros(2,2,2);
Pxxy(1,1,1) = 0.25;
Pxxy(1,2,2) = 0.25;
Pxxy(2,1,2) = 0.25;
Pxxy(2,2,1) = 0.25;

lat = lattice3d();
lat = calc_pe(lat,Pxxy,@Hcs);

fprintf(1,'XOR PED\n')
lat.PI

%% XOR PID
% Griffith & Koch (2014) Fig. 6.5
Pxxy = zeros(2,2,2);
Pxxy(1,1,1) = 0.25;
Pxxy(1,2,2) = 0.25;
Pxxy(2,1,2) = 0.25;
Pxxy(2,2,1) = 0.25;
fprintf(1,'XOR PID\n')
compare_ped(Pxxy);

%% AND PED
Pxxy = zeros(2,2,2);
Pxxy(1,1,1) = 0.25;
Pxxy(1,2,1) = 0.25;
Pxxy(2,1,1) = 0.25;
Pxxy(2,2,2) = 0.25;
lat = lattice3d();
lat = calc_pe(lat,Pxxy,@Hcs);
fprintf(1,'AND PED\n')
lat.PI
%% AND PID
Pxxy = zeros(2,2,2);
Pxxy(1,1,1) = 0.25;
Pxxy(1,2,1) = 0.25;
Pxxy(2,1,1) = 0.25;
Pxxy(2,2,2) = 0.25;
fprintf(1,'AND PID\n')
lat = compare_ped(Pxxy);

%% AND PID - SWITCH OUTPUT
Pxxy = zeros(2,2,2);
Pxxy(1,1,1) = 0.25;
Pxxy(1,2,1) = 0.25;
Pxxy(2,1,1) = 0.25;
Pxxy(2,2,2) = 0.25;
fprintf(1,'AND PID SWITCH OUTPUT\n')
compare_ped(permute(Pxxy,[1 3 2]));

%% SUM PID 
Pxxy = zeros(2,2,3);
Pxxy(1,1,1) = 0.25;
Pxxy(1,2,2) = 0.25;
Pxxy(2,1,2) = 0.25;
Pxxy(2,2,3) = 0.25;
fprintf(1,'SUM PID\n')
compare_ped(Pxxy);

%% W&B FIGURE 4A
Pxxy = zeros(2,2,3);
Pxxy(1,1,1) = 1/3;
Pxxy(1,2,2) = 1/3;
Pxxy(2,1,3) = 1/3;

fprintf(1,'Williams and Beer Fig 4A\n')
compare_ped(Pxxy);

%% W&B FIGURE 4B
Pxxy = zeros(2,2,3);
Pxxy(1,1,1) = 0.25;
Pxxy(1,2,2) = 0.25;
Pxxy(2,2,2) = 0.25;
Pxxy(2,1,3) = 0.25;

fprintf(1,'Williams and Beer Fig 4B\n')
compare_ped(Pxxy);

%% IMPERFECTRDN : Griffith et al. (2014) Fig. 3
Pxxy = zeros(2,2,2);
Pxxy(1,1,1) = 0.4;
Pxxy(1,2,1) = 0.1;
Pxxy(2,2,2) = 0.5;
% Pxxy(1,1,1) = 0.499;
% Pxxy(1,2,1) = 0.001;
% Pxxy(2,2,2) = 0.5;

fprintf(1,'IMPERFECTRDN Griffith (2014) Fig 3\n');
compare_ped(Pxxy);

%% RNDXOR : Griffith et al. (2014) Fig. 2
% Griffith et al. (2012) Fig. 8 
Pxxy = zeros(4,4,4);
p = 1/8;
Pxxy(1,1,1) = p;
Pxxy(1,2,2) = p;
Pxxy(2,1,2) = p;
Pxxy(2,2,1) = p;
Pxxy(3,3,3) = p;
Pxxy(3,4,4) = p;
Pxxy(4,3,4) = p;
Pxxy(4,4,3) = p;

fprintf(1,'RNDXOR Griffith (2014) Fig 2\n');
compare_ped(Pxxy);

%% SUBTLE : Griffith et al. (2014) Fig. 4
Pxxy = zeros(2,2,4);
Pxxy(1,1,1) = 1/3;
Pxxy(1,2,2) = 1/3;
Pxxy(2,2,4) = 1/3;

fprintf(1,'SUBTLE Griffith (2014) Fig 4\n');
compare_ped(Pxxy);


%% RDNUNQXOR : Griffith & Koch (2014) Fig. 6.12
% Out of memory error in maximum entropy function from DIT
% % Pxxy = zeros(8,8,16);
% % p = 1./32;
% % Pxxy(1,1,1) = p;
% % Pxxy(1,2,2) = p;
% % Pxxy(2,1,2) = p;
% % Pxxy(2,2,1) = p;
% % 
% % Pxxy(1,3,3) = p;
% % Pxxy(1,4,4) = p;
% % Pxxy(2,3,4) = p;
% % Pxxy(2,4,3) = p;
% % 
% % Pxxy(3,1,5) = p;
% % Pxxy(3,2,6) = p;
% % Pxxy(4,1,6) = p;
% % Pxxy(4,2,5) = p;
% % 
% % Pxxy(3,3,7) = p;
% % Pxxy(3,4,8) = p;
% % Pxxy(4,3,8) = p;
% % Pxxy(4,4,7) = p;
% % 
% % Pxxy(5,5,9) = p;
% % Pxxy(5,6,10) = p;
% % Pxxy(6,5,10) = p;
% % Pxxy(6,6,9) = p;
% % 
% % Pxxy(5,7,11) = p;
% % Pxxy(5,8,12) = p;
% % Pxxy(6,7,12) = p;
% % Pxxy(6,8,11) = p;
% % 
% % Pxxy(7,5,13) = p;
% % Pxxy(7,6,14) = p;
% % Pxxy(8,5,14) = p;
% % Pxxy(8,6,13) = p;
% % 
% % Pxxy(7,7,15) = p;
% % Pxxy(7,8,16) = p;
% % Pxxy(8,7,16) = p;
% % Pxxy(8,8,15) = p;
% % 
% % fprintf(1,'RDNUNQXOR Griffith & Koch (2014) Fig. 6.12\n');
% % compare_ped(Pxxy);

