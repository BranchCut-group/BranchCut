%% Setup
clear all; close all; clc;

%% Import data
cols = 1959;
rows = 2065;
pha = readBinFile('data/Bam_earthquake/20031203_20040211.cut.pha', cols);

figure(1);
imagesc(pha, [-pi pi]);
colormap(jet);

%% Goldstein Filtering

z = exp(1i*pha); % Generate complex image from phase
z = z./abs(z); % Normalize
zGold = goldstein_filt(z, 50, 0.8);

phaGold = angle(zGold); % Convert back into phase

figure(2);
imagesc(phaGold, [-pi pi]);
colormap(jet);

%% Write Goldstein filtered phase

writeBinFile('data/Bam_earthquake/20031203_20040211.cut.pha.gold', phaGold, 1, 'ieee-be', 'float32')