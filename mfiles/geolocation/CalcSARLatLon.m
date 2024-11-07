function [latP, lonP, hP] = CalcSARLatLon(sample, line, height, sarParameters,  ...
                                          elliSemiMajor, elliFlat)
%
% Creator John Peter Merryman Boncori - sarmap 
% Date: 11 Jan 2017
%
% Usage: 
%
% [latP, lonP, hP] = CalcSARLatLon(sample, line, height, sarParameters, [elliSemiMajor], [elliFlat]);
%
% Compute latitude and longitude for given SAR image coordinates (sample and line).
%
% sample     : (input)  point range column number in SAR image (floating point value)
% line       : (input)  point azimuth line number in SAR image (floating point value)
% height     : (input)  height of point above ellipsoid (m)
% sarParams  : (input)  SLC parameter struct (read in with ParseSMLpar or ParseISPpar)
% ellSemiMaj : (input)  ellipsoid semi major axis (m) (default = WGS84 parameters)
% ellFlat    : (input)  ellipsoid flattening (default = (default = WGS84 parameters))
% latP       : (output) latitude of point on ground (decimal degrees)
% lonP       : (output) longitude of point on ground (decimal degrees)
% hP         : (output) height of point above ellipsoid (same as input height hopefully!)
%
% Example:
%     S = ParseSMLpar('burst_IW2_5_slc.sml');
%     [lat, lon] = CalcSARLatLon(11256.163030669417, 1300.3576101923582, 384.551, S);
%

if nargin < 4
  error('Insufficient number of input parameters, exiting.');
end

if nargin < 5
  WGS84_MAJ  = 6378137; 
  elliSemiMajor = WGS84_MAJ; 
end

if nargin < 7
  WGS84_FLAT = 1/298.257223563;
  elliFlat = WGS84_FLAT; 
end                                         

% Calc azimuth time
zeroDopplerTime = sarParameters.tstart + (line - 1) / sarParameters.PRF;

% dateStart = sec2date(sarParameters.tstart + date2sec(upper(datestr(sarParameters.startTimeDateStr,'dd-mmm-yyyy'))));
% date = sec2date(zeroDopplerTime + date2sec(upper(datestr(sarParameters.startTimeDateStr,'dd-mmm-yyyy'))));
% fprintf('DEBUG: line          = %d\n',  line);
% fprintf('DEBUG: tzd           = %.20f\n',  zeroDopplerTime);
% fprintf('DEBUG: tzd - tstart  = %.20f\n',  zeroDopplerTime - sarParameters.tstart);
% fprintf('DEBUG: taz           = %s\n',  date);
% fprintf('DEBUG: dateStart      = %s\n',  dateStart);

% Calc range
SOL = physconst('LightSpeed');
range = sarParameters.Rnear + (sample - 1) * SOL / 2 / sarParameters.fs;

% Compute satellite position at zeroDoppler time
if ~isfield(sarParameters, 'satelliteTrajectory')
    tstop = sarParameters.tstart + (sarParameters.nli - 1) / sarParameters.PRF;
    [XPoly, VPoly] = CalcSatTrajectory(sarParameters.stateVectors, sarParameters.tstart, tstop);
else
    XPoly = sarParameters.satelliteTrajectory.XPoly;
    VPoly = sarParameters.satelliteTrajectory.VPoly;
end
XSat = polyvalSV(XPoly, zeroDopplerTime); 
VSat = polyvalSV(VPoly, zeroDopplerTime); 

% Compute approximate solution using nominal look angle
if ~isfield(sarParameters, 'sceneCentreLookAngle')
    nominalLookAngle = 30; % nominal look angle in degrees
else
    nominalLookAngle = sarParameters.sceneCentreLookAngle;
end

% Assign look angle sign (positive for right-looking)
if isfield(sarParameters, 'lookDirection')
    if strcmp(sarParameters.lookDirection, 'RIGHT')
        signCoeff = 1;
    else
        signCoeff = -1;        
    end
else
    warning('Look direction not specified, assuming SAR to be right looking.');
    signCoeff = 1;
end
theta0 = deg2rad(nominalLookAngle)*signCoeff;

% LoS in SAT1 coordinate system (C&W p. 553)
losSat = [-cos(theta0) sin(theta0) 0]'; 
% SAT1 basis in ECEF coordinates
Xhat_ecef = XSat / norm(XSat);
Zhat_ecef = VSat / norm(VSat);
Yhat_ecef = cross(Zhat_ecef, Xhat_ecef);
% Transform from SAT1 to ECEF coords
M = [standing(Xhat_ecef) standing(Yhat_ecef) standing(Zhat_ecef)];
Ulos = M*losSat; 

% fprintf('DEBUG: XSat(1) = %15.12e; XSat(2) = %15.12e; XSat(3) = %15.12e\n',  XSat(1), XSat(2), XSat(3));
% fprintf('DEBUG: VSat(1) = %15.12e; VSat(2) = %15.12e; VSat(3) = %15.12e\n',  VSat(1), VSat(2), VSat(3));
% fprintf('DEBUG: Xhat_ecef(1) = %15.12e; Xhat_ecef(2) = %15.12e; Xhat_ecef(3) = %15.12e\n',  Xhat_ecef(1), Xhat_ecef(2), Xhat_ecef(3));
% fprintf('DEBUG: Yhat_ecef(1) = %15.12e; Yhat_ecef(2) = %15.12e; Yhat_ecef(3) = %15.12e\n',  Yhat_ecef(1), Yhat_ecef(2), Yhat_ecef(3));
% fprintf('DEBUG: Zhat_ecef(1) = %15.12e; Zhat_ecef(2) = %15.12e; Zhat_ecef(3) = %15.12e\n',  Zhat_ecef(1), Zhat_ecef(2), Zhat_ecef(3));
% fprintf('DEBUG: Ulos(1) = %15.12e; Ulos(2) = %15.12e; Ulos(3) = %15.12e\n',  Ulos(1), Ulos(2), Ulos(3));


% Find intersection of the line from satellite to Earth centre and a point on the ellipsoid.
% C&W, pp. 558-559
elliSemiMinor = elliSemiMajor*(1 - elliFlat);
elliEpsilon = (elliSemiMajor/elliSemiMinor)^2  - 1; % second eccentricity squared
elliEccSquared = elliFlat*(2-elliFlat);

F    = (dot(XSat,Ulos) + elliEpsilon*XSat(3)*Ulos(3)) / (1 + elliEpsilon*Ulos(3)^2);
G    = (dot(XSat,XSat) - elliSemiMajor^2 + elliEpsilon*XSat(3)^2) / (1 + elliEpsilon*Ulos(3)^2);
R    = -F - sqrt(F^2 - G);
X0   = XSat + R.*Ulos;

% fprintf('DEBUG: R = %15.12e\n',  R);
% 
% fprintf('\nDEBUG: Radar geometry params:\n');
% fprintf('DEBUG: nli    = %15.12e\n',  sarParameters.nli);
% fprintf('DEBUG: nsa    = %15.12e\n',  sarParameters.nsa);
% fprintf('DEBUG: Rnear  = %15.12e\n',  sarParameters.Rnear);
% fprintf('DEBUG: dr     = %15.12e\n',  sarParameters.dr);
% fprintf('DEBUG: theta  = %15.12e deg\n',  nominalLookAngle);
% fprintf('DEBUG: tstart = %s\n',  sarParameters.startTimeDateStr);
% fprintf('DEBUG: PRF    = %15.12e\n',  sarParameters.PRF);

[lat0, lon0, h0] = xyz2ell2(X0(1), X0(2), X0(3), elliSemiMajor, elliEccSquared);
% fprintf('\nDEBUG: Iteration starting point\n');
% fprintf('DEBUG: sample = %.12e ; line  = %.12e ; height = %.12e\n', sample, line, height);
% fprintf('DEBUG: X0(1)  = %.12e ; X0(2) = %.12e ; X0(3)  = %.12e\n', X0(1), X0(2), X0(3));
% fprintf('DEBUG: Lat    = %.12e ; Lon   = %.12e ; height = %.12e;\n', rad2deg(lat0), rad2deg(lon0), h0);

% Solve with Newton method
MAX_ITER = 150; % max. number of iterations
eps = 1e-6;     % convergence epsilon
Xi_1 = [0;0;0]; % inits
Xi = X0;
iter  = 1;

% Convert m to km to avoid ill conditioning
scaleFac = 1e-03;
VSat_km = VSat*scaleFac;
Xi_1_km = Xi_1*scaleFac;
Xi_km = Xi*scaleFac;
XSat_km = XSat*scaleFac;
range_km = range*scaleFac;
a_plus_h_km = (elliSemiMajor + height)*scaleFac;
b_plus_h_km = (elliSemiMinor + height)*scaleFac;

while norm(Xi_km - Xi_1_km) > eps && iter < MAX_ITER

    % Design matrix evaluated at previous solution
    LOS_km = (Xi_km - XSat_km);
    fXi = [dot(VSat_km, LOS_km);...
           dot(LOS_km, LOS_km) - range_km^2;...
           (Xi_km(1)^2 + Xi_km(2)^2) / a_plus_h_km^2 + ...
           (Xi_km(3) / b_plus_h_km)^2 - 1;];

    % Matrix of partial derivatives    
    DfXi = [lying(VSat_km);...
            2*lying(LOS_km);
            2*Xi_km(1)/a_plus_h_km^2 2*Xi_km(2)/a_plus_h_km^2 2*Xi_km(3)/b_plus_h_km^2];      

    % Solve linear system
    DX = DfXi\(-fXi);
%     if iter == 1
%         fprintf('A = \n');
%         for i = 1 : size(DfXi,1)
%             for j = 1 : size(DfXi,2)
%                 fprintf('%.6e ', DfXi(i,j));
%             end
%             fprintf('\n');
%         end
%         fprintf('b = \n');
%         for i = 1 : length(fXi)
%             fprintf('%.6e ', -fXi(i));
%         end
%         fprintf('\n');
%         fprintf('X = \n');
%         for i = 1 : length(DX)
%             fprintf('%.6e ', DX(i));
%         end
%         fprintf('\n');
%         fprintf('rcond() = %.6e\n', rcond(DfXi));
%     end

    % Update solution
    Xi_1_km = Xi_km;
    Xi_km = Xi_1_km + DX;
    
    iter = iter + 1;
end
% Scale result back to meters
Xi = Xi_km / scaleFac;

if iter == MAX_ITER
    latP = nan;
    lonP = nan;
    hP = nan;
    warning('Maximum number of iterations reached without finding a solution.');
else
    [latP, lonP, hP] = xyz2ell2(Xi(1), Xi(2), Xi(3), elliSemiMajor, elliEccSquared);
    latP = rad2deg(latP);
    lonP = rad2deg(lonP);
    if lonP > 180
        lonP = lonP - 360;
    end
end

% Convert llh to salih to check.
% sali = llh2sali([latP lonP hP], sarParameters);
% sampleP = sali(1);
% lineP = sali(2);
% fprintf('DEBUG: \n          Input      Computed     Difference\n');
% fprintf('DEBUG: Sample: %12.6f %12.6f %12.6f\n', sample, sampleP, sample - sampleP);
% fprintf('DEBUG: Line  : %12.6f %12.6f %12.6f\n', line, lineP, line - lineP);
%fprintf('DEBUG: Height: %12.6f %12.6f %12.6f\n', height, hP, height - hP);
  