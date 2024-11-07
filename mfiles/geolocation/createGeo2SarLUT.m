function LUT = createGeo2SarLUT(slcPar, demPar, hgt)
%
% Creator John Peter Merryman Boncori - DTU 
% Date: 28 May 2018
%
% Usage:
%     LUT = createGeo2SarLUT(slcPar, demPar, hgt)
%
% Create a lookup table for SAR image geocoding
%
% slcPar : (input)  Radar geometry parameter file
% demPar : (input)  geocoded gemoetry parameter file
% hgt    : (input)  Raster of ellipsoidal height values in geocoded gemoetry
% LUT    : (output) Complex lookup table in geocoded geometry.
%                   Real part holds range sample numbers.
%                   Imaginary parts the azimuth line numbers.
%

% Compute satellite trajectory (coefficients for orbit interpolation)
tstop = slcPar.tstart + (slcPar.nli - 1) / slcPar.PRF;
[XPoly, VPoly] = CalcSatTrajectory(slcPar.stateVectors, slcPar.tstart, tstop);
slcPar.satelliteTrajectory.XPoly = XPoly;
slcPar.satelliteTrajectory.VPoly = VPoly;

% Put lat,lon.heights in a vector
latVector = demPar.corner_lat - (0 : demPar.nrows - 1) * abs(demPar.post_lat);
lonVector = demPar.corner_lon + (0 : demPar.ncols - 1) * demPar.post_lon;
[LON,LAT] = meshgrid(lonVector,latVector);
llhVector = [LAT(:) LON(:) hgt(:)];

% Convert to radar coordinates
saliVector = llh2sali(llhVector, slcPar);

% Create a complex lookup table (LUT) the size of the output grid.
% The real part will hold the radar range sample numbers.
% The imaginary part will hold the azimuth line numbers.

saLUT = reshape(saliVector(:,1), size(hgt, 1), size(hgt, 2));
liLUT = reshape(saliVector(:,2), size(hgt, 1), size(hgt, 2));
invalidPixels = (saLUT < 1 | saLUT > slcPar.nsa | liLUT < 1 | liLUT > slcPar.nli);
saLUT(invalidPixels) = nan;
liLUT(invalidPixels) = nan;
LUT = complex(saLUT, liLUT);
