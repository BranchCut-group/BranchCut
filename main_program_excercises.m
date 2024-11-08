%% Setup
clear all; close all; clc;
c = 299792458;

%% Load first picture
filepath1 = 'data_2/20191222.slc.par';
binfilepath1 = 'data_2/20191222.slc';

slcPar1 = ParseISPpar(filepath1); % load parameters into a struct called slcPar

col1 = slcPar1.nsa; % slcPar.nsa is used as columns (complex dim)
row1 = slcPar1.nli; % slcPar.li is used as rows

slc1 = readBinFile(binfilepath1,col1,2,'ieee-be','int16');

amplitude1 = abs(slc1);
phase1     = angle(slc1);

%% Load second picture
filepath2 = 'data_2/20191228.slc.par';
binfilepath2 = 'data_2/20191228.slc';

slcPar2 = ParseISPpar(filepath2); % load parameters into a struct called slcPar

col2 = slcPar2.nsa; % slcPar.nsa is used as columns (complex dim)
row2 = slcPar2.nli; % slcPar.li is used as rows

slc2 = readBinFile(binfilepath2,col2,2,'ieee-be','int16');

amplitude2 = abs(slc2);
phase2     = angle(slc2);

%% Calculations
% In order to do a speckle reduction with specific resolution we must
% calculate some of the variables needed.

% Ground speed of the satellite
Vg = CalcVg(slcPar1, slcPar1.nsa/2, slcPar1.nli/2);

% Incidence angle
[look_angle, inc_angle] = CalcLookAndIncidenceAngle(slcPar1, slcPar1.nsa/2, slcPar1.nli/2);

%inc_angle = slcPar1.inc_angle

% Slant-range pixel spacing
r_pix_spc = c/2/slcPar1.fs;
disp(r_pix_spc); % [m]

% Azimuth pixel spacing
az_pix_spc = Vg/slcPar1.PRF;

% Ground range pixel spacing
ground_r_pix_spc = r_pix_spc/sind(inc_angle);

% Calculate looks
target_res = 45; % m
calc_rangeLooks = target_res/ground_r_pix_spc;
calc_azimuthLooks = target_res/az_pix_spc;

disp('target resolution:'),disp(target_res);
disp("range looks:"); disp(calc_rangeLooks);
disp("azimuth looks:"); disp(calc_azimuthLooks);

%% Figures
figure(1)
colormap(gray)
amp_im = imagesc(amplitude1);
xlabel('Slant range (pixels)');
ylabel('Azimuth range (pixels)');

figure(2)
colormap(jet);
phs_im = imagesc(phase1);
xlabel('Slant range');
ylabel('Azimuth range');

%% Form a complex interferogram
% We expect to see a stripey pattern due to the flat surface, with few
% pertubation due to topography and not too much influence from the
% atmosphere
zint = slc1.*conj(slc2);

phaseint = angle(zint);

%% Non-differential interferogram figure
figure(3);
colormap('jet');
zint_corr_im = imagesc(phaseint);
title('Unflattened interferogram')
xlabel('Slant range (pixels)');
ylabel('Azimuth range (pixels)');

%% Remove the flat Earth and topographic phase contributions
% Read flattening phase (.flat file).
flat_phase = readBinFile('data_2/20191222_20191228.flat',col1,1,'ieee-be','float32');

% Apply the phase correction: zint.*exp(-i*flat_phase)
zflat = zint.*exp(-1i*flat_phase);

%% differential interferogram figure
figure(4);
colormap(jet);
geo_amp_im = imagesc(angle(zflat));
xlabel('Slant range (pixels)');
ylabel('Azimuth range (pixels)');

%% Speckle reductions on zflat
azimuthLooks = round(calc_azimuthLooks);
rangeLooks = round(calc_rangeLooks);
h = ones(azimuthLooks, rangeLooks) / azimuthLooks / rangeLooks;

zflatMulti = filter2(h, zflat);
zflatMulti = zflatMulti(1:azimuthLooks:end, 1:rangeLooks:end);

phaseflatMulti = angle(zflatMulti);



%% Goldstein filtering
zflatGoldMulti = goldstein_filt(zflatMulti ./ abs(zflatMulti), 32, 0.5);
zflatGold = goldstein_filt(zflat ./ abs(zflat), 32, 0.5);

phaseflatGold = angle(zflatGold);
phaseflatGoldMulti = angle(zflatGoldMulti);

%% Filtered figures
% Multilooked interferogram figure
figure(5);
subplot(1,3,1);
colormap('jet');
im = imagesc(phaseflatMulti);
xlabel('Slant range (pixels)');
ylabel('Azimuth range (pixels)');
title('Multilook filtering');

% Goldstein filtered interferogram figure
subplot(1,3,2);
title('Goldstein filtering');
colormap('jet');
im = imagesc(phaseflatGold);
xlabel('Slant range (pixels)');
ylabel('Azimuth range (pixels)');
title('Multilook filtering');

% Goldstein and Multilook filtered figure
subplot(1,3,3);
title('Multilooked and Goldstein filtering');
colormap('jet');
im = imagesc(phaseflatGoldMulti);
xlabel('Slant range (pixels)');
ylabel('Azimuth range (pixels)');
title('Multilook filtering');

set(gcf, 'Position',  [100, 250, 1200, 400]);

%% Unwrap interferogram
apix = 100;
rpix = 850;
unwrap_ref = phaseflatGoldMulti(apix,rpix);

residue_map = PhaseResiduesJME(phaseflatGoldMulti);
branchcut_map = BranchCutsJME(residue_map);
phu = FloodFillJME(phaseflatGoldMulti, branchcut_map, [apix, rpix]);

%% Goldstein-Multilook unwrapped interferogram figure
figure(7);
colormap('jet');
unwrap_proc_phase = imagesc(phu);
cb = colorbar;
cb.Label.String = "m";
xlabel('Slant range (pixels)');
ylabel('Azimuth range (pixels)');
hold on;
plot(apix,rpix,'k*');
title('Goldstein-Multilook unwrapped interferogram')

%% Calculate the velocity from the estimated motion

ty = 6/365.25;
sar_lambda = slcPar1.lambda;

del_los = -sar_lambda/(4*pi)*phu./ty;

%% Figure velocity map
figure(8);
colormap('jet');
im = imagesc(del_los);
cb = colorbar;
cb.Label.String = "m/y";
clim([-20 20]);
xlabel('Slant range (pixels)');
ylabel('Azimuth range (pixels)');

%% Geocode the image
demPar = ParseDEMpar('data_2/20191222_gec_lut.dem_par');

LUT = readBinFile('data_2/20191222_gec_lut.bin',demPar.ncols,2);

newLUT = complex((real(LUT) - 1) / rangeLooks + 1, ...
                 (imag(LUT) - 1) / azimuthLooks + 1);

gec_im = geocodeBack(del_los, newLUT);

x = [1:demPar.ncols].*target_res*1e-3;
y = [1:demPar.nrows].*target_res*1e-3;

%% Geocoded map
figure(9);
colormap('jet');
im = imagesc(x,y,gec_im);
cb = colorbar;
cb.Label.String = "m/y";
clim([-20 20]);
xlabel('x dist in km');
ylabel('y dist in km');

%% Calculate Coherence

h = ones(azimuthLooks, rangeLooks) / azimuthLooks / rangeLooks;
num1 = filter2(h, zflat);
den1 = sqrt(filter2(h, abs(slc1).^2));
den2 = sqrt(filter2(h, abs(slc2).^2));
gammaMulti = abs(num1 ./ den1 ./ den2);

gammaMulti = gammaMulti(1:azimuthLooks:end, 1:rangeLooks:end);

%% Coherence map figure
figure(10);
colormap('jet');
im = imagesc(gammaMulti);
xlabel('Slant range (pixels)');
ylabel('Azimuth range (pixels)');
title('Coherence map')
