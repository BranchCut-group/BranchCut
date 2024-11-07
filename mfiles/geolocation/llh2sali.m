function sali = llh2sali(llh, sarParameters, elliParameters)
%
% Creator John Peter Merryman Boncori - DTU 
%         Based on a script from J.J. Mohr (DTU)
% Date: 30 Oct 2017
%
% usage: sali = llh2sali(llh, sarParams, [elliParams])
%
% Compute SAR image coordinates for a given point on ground and
% SAR acquisition geometry.
%
% llh        : (input)  N x 3 vector, each line containing point [lat lon height]
%                       lat, lon in decimal degrees; height in m above reference ellipsoid
% sarParams  : (input)  SLC parameter struct
% ellParams  : (input)  Reference ellipsoid parameters (default = WGS84 parameters)
% sali       : (output) N x 2 vector, each line containing range sample and azimuth line
%                       numbers in SAR image (floating point values, origin in 1,1)
%
% Example:
%     S = ParseISPpar('20070701_20070816.0.isp.slc.par');
%     sali = llh2sali([41.1 12.5 1000], S);
%

if nargin < 2
   error('Insufficient number of input parameters, exiting.');
end

if nargin < 3
   WGS84_MAJ  = 6378137; 
   WGS84_FLAT = 1/298.257223563;
   elliParameters.semiMajorAxis = WGS84_MAJ; 
   elliParameters.flattening = WGS84_FLAT; 
end

% Compute satellite trajectory (State Vector polyonmial coefficients)
tstart = sarParameters.tstart;
tstop = tstart + (sarParameters.nli - 1) / sarParameters.PRF;
if ~isfield(sarParameters, 'satelliteTrajectory')
    [XPoly, VPoly] = CalcSatTrajectory(sarParameters.stateVectors, sarParameters.tstart, tstop);
else
    XPoly = sarParameters.satelliteTrajectory.XPoly;
    VPoly = sarParameters.satelliteTrajectory.VPoly;
end

% Calc point position in ECR xyz
elliEccSquared = elliParameters.flattening * (2 - elliParameters.flattening);
[xPoint, yPoint, zPoint] = ell2xyz(deg2rad(llh(:,1)), deg2rad(llh(:,2)), llh(:,3), ...
    elliParameters.semiMajorAxis, elliEccSquared);

% Find zero-Doppler azimuth time and corresponding slant range
numOfPoints = size(llh, 1);
sali = zeros(numOfPoints, 2);
%parfor (p = 1 : numOfPoints, 8)
for p = 1 : numOfPoints
    XPoint = [xPoint(p); yPoint(p); zPoint(p)];
    sali(p,:) = findZeroDopplerBisect(XPoint, XPoly, VPoly, tstart, tstop);        
end

% Compute pixel numbers (origin = 1,1 following matlab convention)
SOL = physconst('LightSpeed');
rangePixelSpacing = SOL / 2 / sarParameters.fs;
sali(:,1) = 1 + (sali(:,1) - sarParameters.Rnear) / rangePixelSpacing;
sali(:,2) = 1 + (sali(:,2) - sarParameters.tstart) * sarParameters.PRF;

% Convert salih to llh to check.
% p = 1;
% [latP, lonP, heightP] = CalcSARLatLon(sali(p,1), sali(p,2), hP, sarParameters);
% fprintf('DEBUG: \n          Input      Computed     Difference\n');
% fprintf('DEBUG: Lat: %12.6f %12.6f %12.6f\n', llh(p,1), latP, llh(p,1) - latP);
% fprintf('DEBUG: Lon  : %12.6f %12.6f %12.6f\n', llh(p,2), lonP, llh(p,2) - lonP);
% fprintf('DEBUG: Hgt  : %12.6f %12.6f %12.6f\n', llh(p,3), heightP, llh(p,3) - heightP);


function  zeroDopplerRangeTime = findZeroDopplerBisect(pointXYZ, satXYZPoly, satVelPoly, ...
        azimuthStartTime, azimuthStopTime)
%
% Compute zero-Doppler time and slant-range for a given point on ground
%
% pointXYZ             : (input)  ECR Cartesian coordinates of point on ground
% satXYZPoly           : (input)  Polynomial coefficients for satellite position interpolation
% satVelPoly           : (input)  Polynomial coefficients for satellite velocity interpolation
% azimuthStartTime     : (input)  Azimuth start time (seconds)
% zeroDopplerRangeTime : (output) 1 x 2 vector, containing slant-range of closest approach (m)
%                                 and the corresponding zero-Doppler time (seconds) 
%                                  

% Define time interval for zero-Doppler time search
% NOTE: we specify a bit more than the data extent to have a good
%       sensitivity also at the data edges
dataDuration = azimuthStopTime - azimuthStartTime;
searchIntervalStart = azimuthStartTime - dataDuration;
searchIntervalEnd = azimuthStopTime +  dataDuration;

% Init bisection params.
% Start iteration with entire data interval 
searchIntervalDuration = searchIntervalEnd - searchIntervalStart;
aVerySmallTimeInterval = 1e-6; % for convergence criterion
MAX_ITER = 100;                    % max. num. of iterations
iter = 0;   

% While the time search interval is large enough ...
while ((searchIntervalDuration > aVerySmallTimeInterval) && ...
       (iter < MAX_ITER))

    % Bisect time search interval
    searchIntervalDuration = (searchIntervalEnd - searchIntervalStart) / 2;
    trialTime = searchIntervalStart + searchIntervalDuration;

    % Compute satellite position at trial time
    XSat = polyvalSV(satXYZPoly, trialTime);

    % Compute point line of sight vector
    los = pointXYZ - XSat;

    % Compute squint angle
    % NOTE: since we are just interested in the sign of it, we omit
    % normalizing los and vel vectors
    VSat = polyvalSV(satVelPoly, trialTime);
    sinSquintAngle = dot(los, VSat);

    % If the satellite position is beyond the zero Doppler ....
    if (sinSquintAngle < 0)
    % Carry out the next search in the first half of the time interval
        searchIntervalEnd = trialTime;
    else
    % Otherwise search in the second half  
        searchIntervalStart = trialTime;
    end
    
    iter = iter + 1;

end

if iter == MAX_ITER
    error('Maximum number of iterations reached, exiting')
end

% assign zero-Doppler time and slant range at this time
zeroDopplerRangeTime = [norm(los) trialTime];
