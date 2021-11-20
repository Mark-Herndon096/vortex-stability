%% Script to plot and analysus the perturbation evolution of vortex system
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
%%
t_ind = 400;
[max_val, max_ind] = max(s(:,t_ind));
wvs = wvs/.19;

for n = 1:nt
    [max_val2(n),max_ind2(n)] = max(s(:,n));
    wv_max(n) = wvs(max_ind2(n));
end



figure(1)
plot(wvs,log(s(:,t_ind))/tau(t_ind),'k-','LineWidth',1.5)
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
fprintf("Wave number of max A: %6.6f\n",wvs(max_ind)/0.197);

L = zeros(nt,1);
%%
ind2 = 5;
ind3 = 35; 
ind4 = 204;
ind5 = 555;
ind6 = 950;


for k = 1:nt
    %L(k) = log(s(max_ind,k))/tau(k);
    L1(k) = s(max_ind,k);
    L2(k) = s(ind2,k);
    L3(k) = s(ind3,k);
    L4(k) = s(ind4,k);
    L5(k) = s(ind5,k);
    L6(k) = s(ind6,k);
end

leg1 = wvs(max_ind);
leg2 = wvs(ind2);
leg3 = wvs(ind3);
leg4 = wvs(ind4);
leg5 = wvs(ind5);
leg6 = wvs(ind6);

s1 = sprintf('k = %5.4f', leg1);
s2 = sprintf('k = %5.4f', leg2);
s3 = sprintf('k = %5.4f', leg3);
s4 = sprintf('k = %5.4f', leg4);
s5 = sprintf('k = %5.4f', leg5);
s6 = sprintf('k = %5.4f', leg6);

figure(2)
plot(tau,L2,'k-','LineWidth',1.5), hold on
plot(tau,L3,'g-','LineWidth',1.5), hold on
plot(tau,L1,'r-','LineWidth',1.5), hold on
plot(tau,L4,'b-','LineWidth',1.5), hold on
plot(tau,L5,'m-','LineWidth',1.5), hold on
plot(tau,L6,'c-','LineWidth',1.5), hold on
ylim([0.0 7.5])
xlim([0 5])
legend(s2,s3,s1,s4,s5,s6)

figure(3)
plot(tau,wv_max)


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
