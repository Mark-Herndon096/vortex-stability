%% Script to plot the perturbation trajectories of vortex system
clc; clear all; close all;

dir = '/home/markherndon/vortex-stability/DATA/';
fname2 = sprintf('%svortex_trajectories-1000-010.x',dir);
fname = sprintf('%svortex_trajectories-GE-1000-010.x',dir);
fid = fopen(fname,'r','ieee-le');
fid2 = fopen(fname2,'r','ieee-le');
nv = fread(fid,1,'int');
nt = fread(fid,1,'int');

Y   = zeros(nv,nt);
Z   = zeros(nv,nt);
tau = zeros(1,nt);

nv2 = fread(fid2,1,'int');
nt2 = fread(fid2,1,'int');

Y2   = zeros(nv2,nt2);
Z2   = zeros(nv2,nt2);
tau2 = zeros(1,nt2);

for n = 1:nt
	Y(:,n) = fread(fid,nv,'double');
end
for n = 1:nt
	Z(:,n) = fread(fid,nv,'double');
end

tau = fread(fid,nt,'double');

for n = 1:nt
	Y2(:,n) = fread(fid2,nv2,'double');
end
for n = 1:nt
	Z2(:,n) = fread(fid2,nv2,'double');
end

tau2 = fread(fid2,nt2,'double');


%% Plot trajectories
t_ind=2001;

figure(1)

%plot(Y(1,:),Z(1,:),'k-','LineWidth',1.5), hold on

%plot(Y(1,t_ind),Z(1,t_ind),'ko','LineWidth',2), hold on
plot(Y2(2,:),Z2(2,:),'k-','LineWidth',1.5), hold on
plot(Y(2,:),Z(2,:),'r-','LineWidth',1.5), hold on
plot(Y(2,t_ind),Z(2,t_ind),'ro','LineWidth',2), hold on
plot(Y2(2,t_ind),Z2(2,t_ind),'ko','LineWidth',2)
%grid on
ylim([0 10])
xlim([0 2])
legend('Free vortex pair','Ground constraint','Center position','Center position');
%set(L,'Interpreter','LaTeX')