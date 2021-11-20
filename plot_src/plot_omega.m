clc; clear all; close all;


dir = '/home/markherndon/vortex-stability/DATA/';
%fname = sprintf('%sperturbations-1000-001.x',dir);
fname = sprintf('%sOMEGA.x',dir);
fid = fopen(fname,'r','ieee-le');

nk  = fread(fid,1,'int');
wvs = zeros(nk,1);
omg = zeros(nk,1);

a = 0.19;

wvs = fread(fid,nk,'double');
omg = fread(fid,nk,'double');

%%
figure(1)
plot(wvs,omg,'k-','LineWidth',1.5)
ylim([0 1])
xlim([0 3])