% plot as function of predictor dependence.


b = 1/10; 
a = 1/10;

pxy = [(1 + a)/4 (1 - a)/4 (1 - a)/4 (1 + a)/4];
pxz = [(1 + b)/4 (1 - b)/4 (1 - b)/4 (1 + b)/4];

i = 1:21;

pidbroja = zeros(4,length(i));
pidiccs = zeros(4,length(i));
us = zeros(1,length(i));
vs = zeros(1,length(i));
for ii=i
    ii
    u = (6*ii + 94)/800;
    v = 9*(ii - 1)/800;
    
    us(ii) = u;
    vs(ii) = v;
    
    Pxyz = [u pxy(1) - u pxz(1) - u pxy(2) - pxz(1) + u v pxy(3) - v pxz(3) - v pxy(4) - pxz(3) + v];
    Pxyz = reshape(Pxyz,[2 2 2]);
    
    pidbroja(:,ii) = pid_broja(Pxyz);
    lat = calc_pi(lattice2d(),Pxyz,@Iccs_op_me);
    pidiccs(:,ii) = lat.PI;
    
    pz0 = pxz(1) + pxz(2);
    py1 = pxy(2) + pxy(4);
    py0 = 1 - py1;
    c(ii) = (u + v - pz0*(1-py1)) ./ sqrt(py1*(1-py1)*pz0*(1-pz0));
end

ctheory = (3*i - 23)/40;

%%
b = 3/20; 
a = 3/20;

pxy = [(1 + a)/4 (1 - a)/4 (1 - a)/4 (1 + a)/4];
pxz = [(1 + b)/4 (1 - b)/4 (1 - b)/4 (1 + b)/4];

is = 1:21;

pidbroja = zeros(4,length(is));
pidiccs = zeros(4,length(is));
c = zeros(1,length(i));
for ii=1:length(is);
    i = is(ii);
    u = (13*i + 187)/1600;
    v = 17*(i - 1)/1600;
    
    Pxyz = [u pxy(1) - u pxz(1) - u pxy(2) - pxz(1) + u v pxy(3) - v pxz(3) - v pxy(4) - pxz(3) + v];
    Pxyz = reshape(Pxyz,[2 2 2]);
    
    pidbroja(:,ii) = pid_broja(Pxyz);
    lat = calc_pi(lattice2d(),Pxyz,@Iccs_op_me);
    pidiccs(:,ii) = lat.PI; 
    
    pz0 = pxz(1) + pxz(2);
    py1 = pxy(2) + pxy(4);
    py0 = 1 - py1;
    c(ii) = (u + v - pz0*(1-py1)) ./ sqrt(py1*(1-py1)*pz0*(1-pz0));
end


 
%%
figure
subplot(1,2,1)
plot(c,pidbroja)
legend('red','unqX','unqY','syn');
hold on
plot(c,sum(pidbroja,1),'k')

subplot(1,2,2)
plot(c,pidiccs)
hold on
plot(c,sum(pidiccs,1),'k')
%%



l = 0.2;
Nk = 26;
ks = linspace(0,0.25,Nk);

pidbroja = zeros(4,Nk);
pidiccs = zeros(4,Nk);
for ki=2%1:Nk
    P = fit_predpred_binary(ks(ki),l);
    P = P./sum(P(:));
    pidbroja(:,ki) = pid_broja(P);
    lat = calc_pi_wb(lattice2d(),P,@Iccs_pme);
    pidiccs(:,ki) = lat.PI;
end
%%
figure
subplot(1,2,1)
hold all
plot(ks, pidbroja)
plot(ks, sum(pidbroja,1),'k')
legend({'red' 'unqX' 'unqY' 'syn'})
% ylim([-0.2 1])

subplot(1,2,2)
hold all
plot(ks, pidiccs)
plot(ks, sum(pidiccs,1), 'k')
legend({'red' 'unqX' 'unqY' 'syn'})
% ylim([-0.2 1])

suptitle(sprintf('target-pred L = %0.2f',l))
    
   


%% general case
b = 1/10; 
a = 1/10;

pxy = [(1 + a)/4 (1 - a)/4 (1 - a)/4 (1 + a)/4];
pxz = [(1 + b)/4 (1 - b)/4 (1 - b)/4 (1 + b)/4];

minu = max(0, pxz(1)-pxy(2));
minv = max(0, pxz(3)-pxy(4));
maxu = min([ pxy(1) pxz(1) 1+pxz(1)-pxy(2) ]);
maxv = min([ pxy(3) pxz(3) 1+pxz(3)-pxy(4) ]);

minc = 4*(minu + minv) - 1;
maxc = 4*(maxu + maxv) - 1;

Nc = 10;
c = linspace(0,maxc,Nc);
vs = linspace(minv, maxv, Nc);

pidbroja = zeros(4,Nc);
pidiccs = zeros(4,Nc);
ctest = zeros(1,Nc);
us = zeros(1,Nc);
parfor ci=1:Nc
%     ci
    c(ci);
    v = vs(ci);
    u = (c(ci)+1)/4 - v;
    us(ci) = u;
    vs(ci) = v;
    
    Pxyz = [u pxy(1) - u pxz(1) - u pxy(2) - pxz(1) + u v pxy(3) - v pxz(3) - v pxy(4) - pxz(3) + v];
    Pxyz = reshape(Pxyz,[2 2 2]);
   
    pidbroja(:,ci) = pid_broja(Pxyz);
    lat = calc_pi(lattice2d(),Pxyz,@Iccs_op_me);
    pidiccs(:,ci) = lat.PI; 
    
    pz0 = pxz(1) + pxz(2);
    py1 = pxy(2) + pxy(4);
    py0 = 1 - py1;
    ctest(ci) = (u + v - pz0*(1-py1)) ./ sqrt(py1*(1-py1)*pz0*(1-pz0));
end


%%
figure
subplot(1,2,1)
plot(c,pidbroja)
legend('red','unqX','unqY','syn');
hold on
plot(c,sum(pidbroja,1),'k')
grid on

subplot(1,2,2)
plot(c,pidiccs)
hold on
plot(c,sum(pidiccs,1),'k')
grid on



%%
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

%% general case : single parameter
b = 1/10; 
a = 1/10;

pxy = [(1 + a)/4 (1 - a)/4 (1 - a)/4 (1 + a)/4];
pxz = [(1 + b)/4 (1 - b)/4 (1 - b)/4 (1 + b)/4];

xs = {0 9./40};
cranges = {[-0.1 0.1] [1/10 1]};
Ncs = [50 50];
Ncr = length(cranges);

cout = cell(1,Ncr);
pidbrojaout = cell(1,Ncr);
pidiccsout = cell(1,Ncr);

for cri=1:length(cranges)
    Nc = Ncs(cri);
    crange = cranges{cri};
    cs = linspace(crange(1),crange(2),Nc);
    x = xs{cri};
    
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
        ctest(ci) = (u + v - pz0*(1-py1)) ./ sqrt(py1*(1-py1)*pz0*(1-pz0));
    end
    pidbrojaout{cri} = pidbroja;
    pidiccsout{cri} = pidiccs;
    cout{cri} = cs;
end

%%
pidbroja = cell2mat(pidbrojaout);
pidiccs = cell2mat(pidiccsout);
cs = cell2mat(cout);

figure
ax = [];
ax(1) = subplot(1,2,1);
plot(cs,pidbroja(1:3,:))
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

xlim([-0.11 1.01])
ylim([-0.01 0.05])

%%
figure
plot(cout{1},pidiccsout{1})

