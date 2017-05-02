% three variable examples
%% 3 way redundant
Pxxxy = zeros(2,2,2,2);
Pxxxy(1,1,1,1) = 1/2;
Pxxxy(2,2,2,2) = 1/2;

fprintf(1,'RDN: redundant bits\n')
compare(lattice3d(),Pxxxy);

%% XORDUPLICATE : Griffith & Koch (2014) Fig. 6.6
Pxxxy = zeros(2,2,2,2);
p = 1/4;
Pxxxy(1,1,1,1) = p;
Pxxxy(1,2,1,2) = p;
Pxxxy(2,1,2,2) = p;
Pxxxy(2,2,2,1) = p;

fprintf(1,'XORDUPLICATE Griffith & Koch (2014) Fig. 6.6\n')
lat = compare(lattice3d(),Pxxxy);

%% two separate xors
Pxxxy = zeros(2,2,2,4);

p = 1/8;
Pxxxy(1,1,1,1) = p;
Pxxxy(1,1,2,2) = p;
Pxxxy(1,2,1,4) = p;
Pxxxy(1,2,2,3) = p;

Pxxxy(2,1,1,3) = p;
Pxxxy(2,1,2,4) = p;
Pxxxy(2,2,1,2) = p;
Pxxxy(2,2,2,1) = p;

fprintf(1,'Two separate Xors\n')
lat = compare(lattice3d(),Pxxxy);

%% XORLOSES : Griffith & Koch (2014) Fig. 6.7
Pxxxy = zeros(2,2,2,2);
Pxxxy(1,1,1,1) = 1/4;
Pxxxy(1,2,2,2) = 1/4;
Pxxxy(2,1,2,2) = 1/4;
Pxxxy(2,2,1,1) = 1/4;

fprintf(1,'XORLOSSES Griffith & Koch (2014) Fig. 6.7\n')
lat = compare(lattice3d(),Pxxxy);

%% ANDDUPLICATE : Griffith & Koch (2014) Fig. 6.13
Pxxxy = zeros(2,2,2,2);
p = 1/4;
Pxxxy(1,1,1,1) = p;
Pxxxy(1,2,1,1) = p;
Pxxxy(2,1,2,1) = p;
Pxxxy(2,2,2,2) = p;

fprintf(1,'ANDDUPLICATE Griffith & Koch (2014) Fig. 6.13\n')
lat = compare(lattice3d(),Pxxxy);

%% XORUNQ : 2 bit xor, with independent extra bit
Pxxxy = zeros(2,2,2,4);

p = 1/8;
Pxxxy(1,1,1,1) = p;
Pxxxy(1,1,2,2) = p;
Pxxxy(1,2,1,3) = p;
Pxxxy(1,2,2,4) = p;
Pxxxy(2,1,1,3) = p;
Pxxxy(2,1,2,4) = p;
Pxxxy(2,2,1,1) = p;
Pxxxy(2,2,2,2) = p;

fprintf(1,'2 bit XOR with independent extra bit\n')
lat = compare(lattice3d(),Pxxxy);

%% Bertschinger / Rauh. (x1,x2,x1 XOR x2) copied to output
Pxxxy = zeros(2,2,2,8);
Pxxxy(1,1,1,1) = 0.25;
Pxxxy(1,2,2,4) = 0.25;
Pxxxy(2,1,2,6) = 0.25;
Pxxxy(2,2,1,7) = 0.25;

fprintf(1,'Bertschinger / Rauh XOR copied input\n')
lat = compare(lattice3d(),Pxxxy);

%% 3 way EVEN PARITY

Pxxxy = zeros(2,2,2,2);
Pxxxy(1,1,1,1) = 1/8;
Pxxxy(1,1,2,2) = 1/8;
Pxxxy(1,2,1,2) = 1/8;
Pxxxy(1,2,2,1) = 1/8;
Pxxxy(2,1,1,2) = 1/8;
Pxxxy(2,1,2,1) = 1/8;
Pxxxy(2,2,1,1) = 1/8;
Pxxxy(2,2,2,2) = 1/8;

fprintf(1,'3-way EVEN PARITY\n')
lat = compare(lattice3d(),Pxxxy);


%% 3 way AND
Pxxxy = zeros(2,2,2,2);
Pxxxy(1,1,1,1) = 1/8;
Pxxxy(1,1,2,1) = 1/8;
Pxxxy(1,2,1,1) = 1/8;
Pxxxy(1,2,2,1) = 1/8;
Pxxxy(2,1,1,1) = 1/8;
Pxxxy(2,1,2,1) = 1/8;
Pxxxy(2,2,1,1) = 1/8;
Pxxxy(2,2,2,2) = 1/8;

fprintf(1,'3-way AND\n')
lat = compare(lattice3d(),Pxxxy);

%% XORMULTICOAL : Griffith & Koch (2014) Fig. 6.14
Pxxxy = zeros(4,4,4,2);
p = 1/8;
Pxxxy(1,1,1,1) = p;
Pxxxy(4,3,3,1) = p;
Pxxxy(3,4,2,1) = p;
Pxxxy(2,2,4,1) = p;

Pxxxy(3,3,1,2) = p;
Pxxxy(2,1,3,2) = p;
Pxxxy(1,2,2,2) = p;
Pxxxy(4,4,4,2) = p;

fprintf(1,'XORMULTICOAL Griffith & Koch (2014) Fig. 6.14\n')
lat = compare(lattice3d(),Pxxxy);

%% PARITYRDNRDN : Griffith & Koch (2014) (arXiv version) Fig. 13
% Pxxxy = zeros(8,8,8,8);
% p = 1/32;
% 
% Pxxxy(1,1,1,1) = p;
% Pxxxy(1,1,2,2) = p;
% Pxxxy(1,2,1,2) = p;
% Pxxxy(1,2,2,1) = p;
% Pxxxy(2,1,1,2) = p;
% Pxxxy(2,1,2,1) = p;
% Pxxxy(2,2,1,1) = p;
% Pxxxy(2,2,2,2) = p;
% 
% Pxxxy(3,3,3,3) = p;
% Pxxxy(3,3,4,4) = p;
% Pxxxy(3,4,3,4) = p;
% Pxxxy(3,4,4,3) = p;
% Pxxxy(4,3,3,4) = p;
% Pxxxy(4,3,4,3) = p;
% Pxxxy(4,4,3,3) = p;
% Pxxxy(4,4,4,4) = p;
% 
% Pxxxy(5,5,5,5) = p;
% Pxxxy(5,5,6,6) = p;
% Pxxxy(5,6,5,6) = p;
% Pxxxy(5,6,6,5) = p;
% Pxxxy(6,5,5,6) = p;
% Pxxxy(6,5,6,5) = p;
% Pxxxy(6,6,5,5) = p;
% Pxxxy(6,6,6,6) = p;
% 
% Pxxxy(7,7,7,7) = p;
% Pxxxy(7,7,8,8) = p;
% Pxxxy(7,8,7,8) = p;
% Pxxxy(7,8,8,7) = p;
% Pxxxy(8,7,7,8) = p;
% Pxxxy(8,7,8,7) = p;
% Pxxxy(8,8,7,7) = p;
% Pxxxy(8,8,8,8) = p;
% 
% fprintf(1,'PARITYRDNRDN Griffith & Koch (arXiv version) Fig. 13\n')
% lat = compare(lattice3d(),Pxxxy);
