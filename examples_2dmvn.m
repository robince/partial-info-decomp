%% Barrett (2015) Figure 3, left column
a = 0.5;
c = 0.5;
bax = linspace(-0.49,0.99,50);

% % a = 0.5;
% % c = -0.5;
% % bax = linspace(-0.99,0.49,50);
vs = [1 1 1];

Nb = length(bax);

lat = lattice2d();

mmiR = zeros(1,Nb);
mmiS = zeros(1,Nb);
mmiY = zeros(1,Nb);
mmiZ = zeros(1,Nb);
ccsR = zeros(1,Nb);
ccsY = zeros(1,Nb);
ccsZ = zeros(1,Nb);
ccsS = zeros(1,Nb);
for bi=1:Nb
    bi;
    b = bax(bi);
    C = [1 b a; b 1 c; a c 1];
    
    lat = calc_pi_mvn(lat, C, vs, @Iccs_mvn_P2);
    ccsR(bi) = lat.PI(1);
    ccsY(bi) = lat.PI(2);
    ccsZ(bi) = lat.PI(3);
    ccsS(bi) = lat.PI(4);

    lat = calc_pi_mvn(lat, C, vs, @Immi_mvn);
    mmiR(bi) = lat.PI(1);
    mmiY(bi) = lat.PI(2);
    mmiZ(bi) = lat.PI(3);
    mmiS(bi) = lat.PI(4);
end

%%
figure
ax = [];
ax(1) = subplot(2,2,1);
hold all
plot(bax,mmiR)
plot(bax,mmiY)
plot(bax,mmiZ)
plot(bax,mmiS)
plot(bax,mmiR+mmiZ+mmiY+mmiS,'k')
legend('Red','Unq1','Unq2','Syn')
title('Immi')
grid on

ax(2) = subplot(2,2,2);
hold all
plot(bax,ccsR)
plot(bax,ccsY)
plot(bax,ccsZ)
plot(bax,ccsS)
plot(bax,ccsR+ccsZ+ccsY+ccsS,'k')
legend('Red','Unq1','Unq2','Syn')
title('Iccs')
grid on

linkaxes(ax,'xy')
ylim([-0.2 1])

%% Barrett (2015) Figure 3, right column
a = 0.4;
c = 0.6;
bax = linspace(-0.48,0.94,50);

% a = 0.25;
% c = 0.75;
% bax = linspace(-0.45,0.82,50);

% a = 0.25;
% c = -0.75;
% bax = linspace(-0.82,0.45,50);

vs = [1 1 1];

Nb = length(bax);

lat = lattice2d();
ccsR = zeros(1,Nb);
mmiR = zeros(1,Nb);
ccsS = zeros(1,Nb);
mmiS = zeros(1,Nb);
ccsY = zeros(1,Nb);
mmiY = zeros(1,Nb);
ccsZ = zeros(1,Nb);
mmiZ = zeros(1,Nb);
for bi=1:Nb
    bi;
    b = bax(bi);
    C = [1 b a; b 1 c; a c 1];
    
    lat = calc_pi_mvn(lat, C, vs, @Iccs_mvn_P2);
    ccsR(bi) = lat.PI(1);
    ccsY(bi) = lat.PI(2);
    ccsZ(bi) = lat.PI(3);
    ccsS(bi) = lat.PI(4);

    lat = calc_pi_mvn(lat, C, vs, @Immi_mvn);
    mmiR(bi) = lat.PI(1);
    mmiY(bi) = lat.PI(2);
    mmiZ(bi) = lat.PI(3);
    mmiS(bi) = lat.PI(4);
end

%%
ax = [];
ax(1) = subplot(2,2,3);
cla
hold all
plot(bax,mmiR)
plot(bax,mmiY)
plot(bax,mmiZ)
plot(bax,mmiS)
plot(bax,mmiR+mmiZ+mmiY+mmiS,'k')
legend('Red','Unq1','Unq2','Syn')
title('Immi')
grid on

ax(2) = subplot(2,2,4);
cla
hold all
plot(bax,ccsR)
plot(bax,ccsY)
plot(bax,ccsZ)
plot(bax,ccsS)
plot(bax,ccsR+ccsZ+ccsY+ccsS,'k')
legend('Red','Unq1','Unq2','Syn')
title('Iccs')
grid on

linkaxes(ax,'xy')
ylim([-0.2 1])