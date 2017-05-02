%% general case : single parameter
b = 1/10; 
a = 1/10;

pxy = [(1 + a)/4 (1 - a)/4 (1 - a)/4 (1 + a)/4];
pxz = [(1 + b)/4 (1 - b)/4 (1 - b)/4 (1 + b)/4];

% xs = {0 9./40};
% cranges = {[-0.1 0.1] [1/10 1]};
% Ncs = [50 50];
% Ncr = length(cranges);

crange = [-0.1 0.1];
% fix x=0 for this range
x = 0;

Nc = 50;
cs = linspace(crange(1),crange(2),Nc);

pidbroja = zeros(4,Nc);
pidiccs = zeros(4,Nc);
parfor ci=1:Nc
    c = cs(ci);
    pyz = [(1 + c)/4 (1 - c)/4 (1 - c)/4 (1 + c)/4];
    
    Pxyz = [pyz(1)-x ...
        pxz(2) - pxy(2) + pxz(1) - pyz(1) + x ...
        pxz(1) - pyz(1) + x ...
        pxy(2) - pxz(1) + pyz(1) - x ...
        x ...
        pxy(3) - x ...
        pxz(3) - x ...
        pxy(4) - pxz(3) + x ];
    Pxyz = reshape(Pxyz,[2 2 2]);
    
    pidbroja(:,ci) = pid_broja(Pxyz);
    lat = calc_pi(lattice2d(),Pxyz,@Iccs_op_me);
    pidiccs(:,ci) = lat.PI;
    
    pz0 = pxz(1) + pxz(2);
    py1 = pxy(2) + pxy(4);
    py0 = 1 - py1;
end

%%
figure
ax = [];
ax(1) = subplot(1,2,1);
pidm = mean(pidbroja(1:3,:),2);

plot(cs,pidm(:,ones(1,Nc)))
legend('red','unqX','unqY');
hold on
% plot(cs,sum(pidbroja,1),'k')
grid on

ax(2) = subplot(1,2,2);
plot(cs,pidiccs(1:3,:))
hold on
% plot(cs,sum(pidiccs,1),'k')
grid on

linkaxes(ax,'xy')



