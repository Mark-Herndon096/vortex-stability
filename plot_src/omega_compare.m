clc; clear all; close all;


nk = 400;
ka = linspace(0.001, 1.5, nk);

omega_1 = zeros(nk,1);
omega_2 = zeros(nk,1);

for n = 1:nk
	omega_1(n) = cut_off(ka(n));
    omega_s(n) = omega_1(n)*ka(n)^2;
	omega_2(n) = self_induct(ka(n));
end

%% Plot omegas
figure(1)
plot(ka,omega_s,'k-','LineWidth',1.5), hold on
plot(ka,omega_2,'r-','LineWidth',1.5)
ylim([0 1])
xlim([0 3])
%%
% figure(2)
% semilogx(ka,omega_1,'k-','LineWidth',1.5)
% xlim([0.1 10])
% grid on

%% Section for eigenvalues from different version of omega
alpha_crow = zeros(nk,1);
alpha_saf = zeros(nk,1);
a = 0.098;
for n = 1:nk
   psi = psi_func(ka(n)/a);
   chi = chi_func(ka(n)/a);
   omega = omega_2(n)/a^2;
   t1 = 1.0 - psi + omega;
   t2 = 1.0 + chi - omega;
   alpha_saf(n,1) = sqrt(t1*t2);
end
alpha_saf(:,1) = real(alpha_saf(:,1));
a = 0.098;
for n = 1:nk
   psi = psi_func(ka(n)/a);
   chi = chi_func(ka(n)/a);
   omega = (ka(n)^2)*omega_1(n)/a^2;
   t1 = 1.0 - psi + omega;
   t2 = 1.0 + chi - omega;
   alpha_crow(n,1) = sqrt(t1*t2);
end
alpha_crow(:,1) = real(alpha_crow(:,1));
a1 = 0.098;
a2 = 0.098;
%%
figure(3)
plot(ka(3:end)/a1,alpha_saf(3:end),'ko','LineWidth',1.5), hold on
plot(ka/a2,alpha_crow,'r-','LineWidth',1.5)
xlim([0.0 4])
ylim([0 1])
legend('Dispersion relation','Cut off method');
%set(L1,'Interpreter','LaTeX')
%%
[max_val, max_ind] = max(alpha_saf(:,1));
[max_val2, max_ind2] = max(alpha_crow(:,1));

fprintf('Max alpha = %6.5f for saffman rate\n',max_val)
fprintf('Max alpha = %6.5f for crow    rate\n',max_val2)
fprintf('Max beta  = %6.5f for saffman rate\n',ka(max_ind)/a1)
fprintf('Max beta  = %6.5f for Crow rate\n',ka(max_ind2)/a2)


%% Read in omega from simulations
dir = '/home/markherndon/vortex-stability/DATA/';
%fname = sprintf('%sperturbations-1000-001.x',dir);
fname = sprintf('%sOMEGA.x',dir);
fid = fopen(fname,'r','ieee-le');

nk2  = fread(fid,1,'int');
wvs2 = zeros(nk2,1);
omg2 = zeros(nk2,1);

a = 0.19;

wvs2 = fread(fid,nk2,'double');
omg2 = fread(fid,nk2,'double');

%%
figure(4)
plot(wvs2,omg2,'bo','LineWidth',1.5), hold on
plot(ka,omega_s,'k-','LineWidth',1.5), hold on
plot(ka,omega_2,'r-','LineWidth',1.5)
ylim([0 0.52])
xlim([0 1.5])


