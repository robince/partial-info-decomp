
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
ccsR = zeros(1,Nb);
mmiR = zeros(1,Nb);
ccsS = zeros(1,Nb);
mmiS = zeros(1,Nb);
ccsY = zeros(1,Nb);
mmiY = zeros(1,Nb);
ccsZ = zeros(1,Nb);
mmiZ = zeros(1,Nb);
for bi=1:Nb
    bi
    b = bax(bi);
    C = [1 b a; b 1 c; a c 1];
    
    lat = calc_pi_mvn(lat, C, vs, @Iccs_mvn);
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
subplot(1,2,1)
plot(bax,ccsS,'r')
hold all
% plot(bax,ccsR,'r--')
plot(bax, mean(ccsR)*ones(size(bax)),'r--')
plot(bax,mmiS,'b')
plot(bax,mmiR,'b--')
xlim([-1 1])
hline(0,':k')
vline(-0.5,':k')
ylim([-0.2 1])
[mean(ccsY) mmiY(1)]
[mean(ccsZ) mmiZ(1)]
legend('ccsS','ccsR','mmiS','mmiR')
%% Barrett (2015) Figure 3, right column
a = 0.25;
c = 0.75;
bax = linspace(-0.45,0.82,50);

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
    bi
    b = bax(bi);
    C = [1 b a; b 1 c; a c 1];
    
    lat = calc_pi_mvn(lat, C, vs, @Iccs_mvn);
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
% figure
subplot(1,2,2)
cla
plot(bax,ccsS,'r')
hold all
% plot(bax,ccsR,'r--')
plot(bax, mean(ccsR)*ones(size(bax)),'r--')
plot(bax,mmiS,'b')
plot(bax,mmiR,'b--')
xlim([-1 1])
ylim([-0.2 1])
hline(0,':k')
vline(-0.4529,':k')
vline(0.8279,':k')
[mean(ccsY) mmiY(1)]
[mean(ccsZ) mmiZ(1)]
% legend('ccsS','ccsR','mmiS','mmiR')