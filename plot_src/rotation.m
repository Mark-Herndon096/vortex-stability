%% Self induction function comparison
clc; clear all; close all;

np = 250;

omega      = zeros(1,np);
omega_crow = zeros(1,np);

start_val = 1e-5;
end_val   = 20;

k = linspace(start_val, end_val, np);

b = 1.0; 
a = 0.63;
d = 0.321*a;
eps = 0.642*a;

ka = k*a;
keps = k*eps;
kd = k*d;

for n = 1:np
   omega(n) = self_induct(ka(n))/a^2;
   omega_crow(n) = (k(n)^2)*cut_off(kd(n));
end

%% Plot self induced rotation rates
close all;
figure(1)
plot(k,omega,'k-','LineWidth',2), hold on
plot(k,omega_crow,'k--','LineWidth',2)
grid on
xlim([0 2])
ylim([0 0.35])
