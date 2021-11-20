%% Script to plot optimal perturbation analysis results
clc; clear all; close all;

dir = '/home/markherndon/vortex-stability/DATA/';
%fname = sprintf('%sperturbations-1000-001.x',dir);
fname = sprintf('%sperturbations-1000-005.x',dir);
fid = fopen(fname,'r','ieee-le');

nv  = fread(fid,1,'int');
nt  = fread(fid,1,'int');
nk  = fread(fid,1,'int');
m   = nv*2;
PHI = zeros(m,m,nt,nk);
tau = zeros(1,nt);
s   = zeros(nk,nt);
wvs = zeros(1,nk);
V = zeros(m,nt,nk);

for k = 1:nk
    for n = 1:nt
        for j = 1:m
            PHI(:,j,n,k) = fread(fid,m,'double');
        end
    end
end

tau = fread(fid,nt,'double');
ind = 1;

for n = 1:nt
    s(:,n) = fread(fid,nk,'double');
end


wvs = fread(fid,nk,'double');

for k = 1:nk
    for n = 1:nt
        V(:,n,k) = fread(fid,m,'double');
    end
end
%%
t_ind = 2000;
a = 0.19;
[max_val, max_ind] = max(s(:,t_ind));



figure(1)
plot(wvs/a,log(s(:,t_ind))/tau(t_ind),'k-','LineWidth',1.5)
xlim([0 3])
ylim([0 1])
grid on
yy3 = ylabel('\(\alpha\)');
set(yy3,'Interpreter','LaTeX','FontSize',16);
xx4 = xlabel('\(\beta\)');
set(xx4,'Interpreter','LaTeX','FontSize',16);
tt4 = title('Maximum amplification SVD approach');
set(tt4,'Interpreter','LaTeX','FontSize',14);
fprintf("-----------------------------\n");
fprintf("Amplification factor: %6.6f\n",log(max_val)/tau(t_ind));
fprintf("Wave number of max A: %6.6f\n",wvs(max_ind)/a);

L = zeros(nt,1);

for k = 1:nt
    L(k) = log(s(max_ind,k))/tau(k);
end

figure(2)
plot(tau,L,'r-','LineWidth',1.5)
grid on
ylim([0.05 0.9])

%%

% figure(3)
% for tt = 1:nt
% plot(wvs/0.19,log(s(:,tt))/tau(tt),'k-','LineWidth',1.5)
% xlim([0 3])
% ylim([0 1])
% grid on
% %fprintf("-----------------------------\n");
% %fprintf("Amplification factor: %6.6f\n",log(max_val)/tau(t_ind));
% %fprintf("Wave number of max A: %6.6f\n",wvs(max_ind)/0.197);
% m = 0.05;
% pause(m)
% clf;
% end
