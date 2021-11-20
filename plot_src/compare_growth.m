%% Main program to compare singular values from GE case and non-GE case
clc; clear all; close all;

dir = '/home/markherndon/vortex-stability/DATA/';
a = 0.098;

%% NON GE variable allocation
fname  = sprintf('%sperturbations_2-1000-020.x',dir);
fid    = fopen(fname,'r','ieee-le');

nv  = fread(fid,1,'int');
nt  = fread(fid,1,'int');
nk  = fread(fid,1,'int');
m   = nv*2;
PHI = zeros(m,m,nt,nk);
tau = zeros(1,nt);
s   = zeros(nk,nt);
wvs = zeros(1,nk);
V = zeros(m,nt,nk);
%% Variable allocation for GE case

fname2  = sprintf('%sperturbations_2-GE-1000-010.x',dir);
fid2    = fopen(fname2,'r','ieee-le');

nv2  = fread(fid2,1,'int');
nt2  = fread(fid2,1,'int');
nk2  = fread(fid2,1,'int');
m2   = nv2*2;
PHI2 = zeros(m2,m2,nt2,nk2);
tau2 = zeros(1,nt2);
s2   = zeros(nk2,nt2);
wvs2 = zeros(1,nk2);
V2 = zeros(m,nt,nk);

%fprintf("All variables allocated for data read\n")
%% Read data into variables for non-GE case
% for k = 1:nk
%     for n = 1:nt
%         for j = 1:m
%             PHI(:,j,n,k) = fread(fid,m,'double');
%         end
%     end
% end

%fprintf("Phi matrix data read\n")

tau = fread(fid,nt,'double');

for n = 1:nt
    s(:,n) = fread(fid,nk,'double');
end

wvs = fread(fid,nk,'double')/a;
% for k = 1:nk
%     for n = 1:nt
%         V(:,n,k) = fread(fid,m,'double');
%     end
% end
%% Read data into variables for GE case
% for k = 1:nk2
%     for n = 1:nt2
%         for j = 1:m2
%             PHI2(:,j,n,k) = fread(fid2,m2,'double');
%         end
%     end
% end
%fprintf("Phi2 matrix read\n")
tau2 = fread(fid2,nt2,'double');

for n = 1:nt2
    s2(:,n) = fread(fid2,nk2,'double');
end

wvs2 = fread(fid2,nk2,'double')/a;
% for k = 1:nk2
%     for n = 1:nt2
%         V2(:,n,k) = fread(fid2,m2,'double');
%     end
% end



%% Figures for comparison
t_ind = 24000;

[waves, growth, np] = analytical_crow;
%for n = 1:np
%    growth(n) = exp(growth(n)*tau(t_ind));
%end


%%
 figure(6)
 plot(wvs(1:end),s(1:end,t_ind)/tau(t_ind),'k-','LineWidth',1.5), hold on
 plot(waves,growth/tau(t_ind),'r-','LineWidth',1.5)
 xlim([0.0 1.5])
 legend('SVD approach','Crow''''s theory')
alp = max(log(s(1:end,t_ind))/tau(t_ind));
%%

figure(7)
plot(wvs(1:5:end),log(s(1:5:end,t_ind))/tau(t_ind),'ko','LineWidth',1.5), hold on
plot(waves,growth,'r-','LineWidth',1.5)
%plot(wvs2(1:end),log(s2(1:end,t_ind))/tau(t_ind),'r-','LineWidth',1.5)
xlim([0.0 4])
ylim([0 1])
legend('SVD approach','Crow''''s theory')
alp = max(log(s(1:end,t_ind))/tau(t_ind));











%%
t_ind = 2001;
[sig_max] = max(s(:,t_ind));

figure(1)
plot(wvs(2:end),s(2:end,t_ind)/sig_max,'k-','LineWidth',1.5), hold on
plot(wvs2(2:end),s2(2:end,t_ind)/sig_max,'r-','LineWidth',1.5)
xlim([0.03 1.5])
ylim([0 1.1])
legend('Free vortex pair','Ground constraint')

[max_val, max_ind] = max(s(:,t_ind));
[max_val2, max_ind2] = max(s2(:,t_ind));
k_ind1= max_ind;
k_ind2 = max_ind2;


fprintf("Wave number of maximum growth = %6.5f\n",wvs(k_ind1))
%fprintf("Wave number of maximum growth = %6.5f\n",wvs2(k_ind2))


%% Calculate max wave number as function of time
max_wv = zeros(nt,2);
for n = 1:nt
    [max_val, k_max] = max(s(:,n));
    [max_val2, k_max2] = max(s2(:,n));
    max_wv(n,1) = wvs(k_max);
    max_wv(n,2) = wvs2(k_max2);
end
%%
figure(2)
plot(tau,s(k_ind1,:)/1000,'k-','LineWidth',1.5), hold on
plot(tau2,s2(k_ind2,:)/1000,'r-','LineWidth',1.5)
ylim([0 4])
xlim([0 tau(end)])
legend('Free vortex pair','Ground constraint','location','northwest')
%%

%% Construct max solution
% fprintf("Constructing max solution")
% x0 = V(:,t_ind,k_ind1);
% x02 = V2(:,t_ind,k_ind2);
% q1 = zeros(m,nt);
% q2 = zeros(m2,nt2);
% 
% for n = 1:nt
%     q1(:,n) = PHI(:,:,n,k_ind1)*x0;
% end
% for n = 1:nt2
%     q2(:,n) = PHI2(:,:,n,k_ind2)*x02;
% end
% 
% sq1 = zeros(nt,1);
% sq2 = zeros(nt2,1);
% den1 =  sqrt(q1(1,1)^2 + q1(3,1)^2);
% den2 =  sqrt(q2(1,1)^2 + q2(3,1)^2);
% for n = 1:nt
%     sq1(n) = sqrt(q1(1,n)^2 + q1(3,n)^2)/den1;
% end
% 
% for n = 1:nt2
%     sq2(n) = sqrt(q2(1,n)^2 + q2(3,n)^2)/den2;
% end

% 
% %%
% figure(3)
% plot(tau,max_wv(:,1),'k-','LineWidth',1.5), hold on
% plot(tau,max_wv(:,2),'r-','LineWidth',1.5)
% tau(t_ind)
% 
%  figure(4)
% plot(tau,sq1/10,'k-','LineWidth',1.5), hold on
% plot(tau2,sq2/10,'r-','LineWidth',1.5)
% ylim([0 1.5])
% xlim([0 1.5])




