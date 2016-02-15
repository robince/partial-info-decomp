% two variable examples

%% W&B FIGURE 4A
Pxxy = zeros(2,2,3);
Pxxy(1,1,1) = 1/3;
Pxxy(1,2,2) = 1/3;
Pxxy(2,1,3) = 1/3;

fprintf(1,'Williams and Beer Fig 4A\n')
lat = compare(lattice2d(),Pxxy);

%% W&B FIGURE 4B
Pxxy = zeros(2,2,3);
Pxxy(1,1,1) = 0.25;
Pxxy(1,2,2) = 0.25;
Pxxy(2,2,2) = 0.25;
Pxxy(2,1,3) = 0.25;

fprintf(1,'Williams and Beer Fig 4B\n')
compare(lattice2d(),Pxxy);

%% W&B FIGURE 4 MODIFIED C
Pxxy = zeros(2,2,3);
p = 1/6;
Pxxy(1,1,1) = p;
Pxxy(1,2,2) = p;
Pxxy(2,2,2) = p;
Pxxy(2,1,3) = p;
Pxxy(2,2,1) = p;
Pxxy(1,2,3) = p;

fprintf(1,'Williams and Beer Fig 4 Modified 1\n')
compare(lattice2d(),Pxxy);

%% W&B FIGURE 4 MODIFIED D
Pxxy = zeros(2,2,3);
p = 1/5;
Pxxy(1,1,1) = p;
Pxxy(1,2,2) = p;
Pxxy(2,2,2) = p;
Pxxy(2,1,3) = p;
Pxxy(2,2,1) = p;
% Pxxy(1,2,3) = p;

fprintf(1,'Williams and Beer Fig 4 Modified 2\n')
compare(lattice2d(),Pxxy);

%% one bit copy
Pxxy = zeros(2,2,2);
Pxxy(1,1,1) = 0.25;
Pxxy(1,2,1) = 0.25;
Pxxy(2,1,2) = 0.25;
Pxxy(2,2,2) = 0.25;

fprintf(1,'1 bit copy\n')
compare(lattice2d(),Pxxy);

%% two bit copy (independent bits)
% UNQ : Griffith et al. (2014) Fig. 1
% UNQ : Griffith et al. (2012) Fig. 3
Pxxy = zeros(2,2,4);
Pxxy(1,1,1) = 0.25;
Pxxy(1,2,2) = 0.25;
Pxxy(2,1,3) = 0.25;
Pxxy(2,2,4) = 0.25;

fprintf(1,'2 bit copy (independent)\n')
compare(lattice2d(),Pxxy);

%% 2 bit copy (redundant bits)
Pxxy = zeros(2,2,4);
Pxxy(1,1,1) = 0.25;
Pxxy(1,1,2) = 0.25;
Pxxy(2,2,3) = 0.25;
Pxxy(2,2,4) = 0.25;
fprintf(1,'2 bit copy (redundant)\n')
compare(lattice2d(),Pxxy);

%% Unequal independent copy
Pxxy = zeros(4,2,8);
p = 1/8;
Pxxy(1,1,1) = p;
Pxxy(2,1,2) = p;
Pxxy(3,1,3) = p;
Pxxy(4,1,4) = p;
Pxxy(1,2,5) = p;
Pxxy(2,2,6) = p;
Pxxy(3,2,7) = p;
Pxxy(4,2,8) = p;

fprintf(1,'Unequal independent copy\n')
compare(lattice2d(),Pxxy);

%% OR
Pxxy = zeros(2,2,2);
Pxxy(1,1,1) = 0.25;
Pxxy(1,2,2) = 0.25;
Pxxy(2,1,2) = 0.25;
Pxxy(2,2,2) = 0.25;
fprintf(1,'OR\n')
lat = compare(lattice2d(),Pxxy);

% lat = lattice2d();
% lat = calc_pi_ri(lat,Pxxy,@Iccs);
% latmin = calc_pi_wb(lat,Pxxy,@Imin);
%% XOR
% Griffith et al. (2012) Fig 4
Pxxy = zeros(2,2,2);
Pxxy(1,1,1) = 0.25;
Pxxy(1,2,2) = 0.25;
Pxxy(2,1,2) = 0.25;
Pxxy(2,2,1) = 0.25;
fprintf(1,'XOR\n')
compare(lattice2d(),Pxxy);

%% AND
Pxxy = zeros(2,2,2);
Pxxy(1,1,1) = 0.25;
Pxxy(1,2,1) = 0.25;
Pxxy(2,1,1) = 0.25;
Pxxy(2,2,2) = 0.25;
fprintf(1,'AND\n')
lat = compare(lattice2d(),Pxxy);

%% SUM
Pxxy = zeros(2,2,3);
Pxxy(1,1,1) = 0.25;
Pxxy(1,2,2) = 0.25;
Pxxy(2,1,2) = 0.25;
Pxxy(2,2,3) = 0.25;
fprintf(1,'SUM\n')
lat = compare(lattice2d(),Pxxy);

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
lat = compare(lattice2d(),Pxxy);

%% IMPERFECTRDN : Griffith et al. (2014) Fig. 3
Pxxy = zeros(2,2,2);
Pxxy(1,1,1) = 0.4;
Pxxy(1,2,1) = 0.1;
Pxxy(2,2,2) = 0.5;
% Pxxy(1,1,1) = 0.499;
% Pxxy(1,2,1) = 0.001;
% Pxxy(2,2,2) = 0.5;

fprintf(1,'IMPERFECTRDN Griffith (2014) Fig 3\n');
lat = compare(lattice2d(),Pxxy);
% NB: shows higher {1}{2} than Imin! This is because {2} has unique
% misinformation

%% SUBTLE : Griffith et al. (2014) Fig. 4
Pxxy = zeros(2,2,4);
Pxxy(1,1,1) = 1/3;
Pxxy(1,2,2) = 1/3;
Pxxy(2,2,4) = 1/3;

fprintf(1,'SUBTLE Griffith (2014) Fig 4\n');
compare(lattice2d(),Pxxy);

%% RDN : Griffith et al. (2012) Fig. 2
Pxxy = zeros(2,2,2);
Pxxy(1,1,1) = 0.5;
Pxxy(2,2,2) = 0.5;

fprintf(1,'RDN Griffith (2012) Fig 2\n');
compare(lattice2d(),Pxxy);

%% XORAND : Bertschinger 2014 Table 1
Pxxy = zeros(2,2,3);
Pxxy(1,1,1) = 0.25;
Pxxy(1,2,3) = 0.25;
Pxxy(2,1,3) = 0.25;
Pxxy(2,2,2) = 0.25;

fprintf(1,'XORAND Bertschinger (2014) Table 1\n');
lat = compare(lattice2d(),Pxxy);


%% RDNUNQXOR : Griffith et al. (2012) Fig 12
Pxxy = zeros(8,8,16);
p = 1./32;
Pxxy(1,1,1) = p;
Pxxy(1,2,2) = p;
Pxxy(2,1,2) = p;
Pxxy(2,2,1) = p;

Pxxy(1,3,3) = p;
Pxxy(1,4,4) = p;
Pxxy(2,3,4) = p;
Pxxy(2,4,3) = p;

Pxxy(3,1,5) = p;
Pxxy(3,2,6) = p;
Pxxy(4,1,6) = p;
Pxxy(4,2,5) = p;

Pxxy(3,3,7) = p;
Pxxy(3,4,8) = p;
Pxxy(4,3,8) = p;
Pxxy(4,4,7) = p;



Pxxy(5,5,9) = p;
Pxxy(5,6,10) = p;
Pxxy(6,5,10) = p;
Pxxy(6,6,9) = p;

Pxxy(5,7,11) = p;
Pxxy(5,8,12) = p;
Pxxy(6,7,12) = p;
Pxxy(6,8,11) = p;

Pxxy(7,5,13) = p;
Pxxy(7,6,14) = p;
Pxxy(8,5,14) = p;
Pxxy(8,6,13) = p;

Pxxy(7,7,15) = p;
Pxxy(7,8,16) = p;
Pxxy(8,7,16) = p;
Pxxy(8,8,15) = p;

fprintf(1,'RDNUNQXOR Griffith (2012) Fig 2\n');
compare(lattice2d(),Pxxy);
