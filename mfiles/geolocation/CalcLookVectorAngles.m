function [theta, phi] = CalcLookVectorAngles(sarParameters, sample, line, height, ...
                                                      elliSemiMajor, elliFlat)
%
% Creator John Peter Merryman Boncori - sarmap 
% Date: 11 Jan 2017
%
% Usage: 
%
% [theta, phi] = CalcLookVectorAngles(sarParameters, sample, line, [height], ...
%                           [elliSemiMajor], [elliFlat]);
%
% Compute line-of-sight theta and phi angles for a given pixel. These
% angles follow the GAMMA software convention:
% theta: SAR look-vector elevation angle (PI/2 -> up  -PI/2 -> down)
% phi  : SAR look-vector orientation angle (lv_phi: 0 -> East  PI/2 -> North) 
%
% sarParameters     : (input)  Raw (or SLC) data parameter struct
% sample            : (input)  Raw or SLC range sample number
% line              : (input)  Raw or SLC azimuth line number
% height            : (input)  pixel height (default = 500 m)
% ellSemiMaj        : (input)  ellipsoid semi major axis (m) (default = WGS84 parameters)
% ellFlat           : (input)  ellipsoid flattening (default = (default = WGS84 parameters))
% theta             : (output) Elevation angle (rad)
% phi               : (output) Orientation angle (rad)
%
% Example:
%     S = ParseISPpar('19960410.rslc.par');
%     [theta, phi] = CalcLookVectorAngles(S, S.nsa/2, S.nli/2);
%                                                        
                                                  
if nargin < 3
  error('Insufficient number of input parameters, exiting.');
end

if nargin < 4
  height  = 500; 
end

if nargin < 5
  WGS84_MAJ  = 6378137; 
  elliSemiMajor = WGS84_MAJ; 
end

if nargin < 6
  WGS84_FLAT = 1/298.257223563;
  elliFlat = WGS84_FLAT; 
end   

% Compute coordinates of the point on ground for which the given sample/line
% represent the zero-Doppler radar coordinates
elliEccSquared = elliFlat*(2-elliFlat);
[latP, lonP] = CalcSARLatLon(sample, line, height, sarParameters);

% Convert to Cartesian coordinates.
[xP, yP, zP] = ell2xyz(deg2rad(latP), deg2rad(lonP), height, elliSemiMajor, elliEccSquared);
PointXYZ = [xP; yP; zP];

% Compute satellite position and velocity polynomials from state vectors
azimuthTime = sarParameters.tstart + (line - 1) / sarParameters.PRF;
if ~isfield(sarParameters, 'satelliteTrajectory')
    tstop = sarParameters.tstart + (sarParameters.nli - 1) / sarParameters.PRF;
    XPoly = CalcSatTrajectory(sarParameters.stateVectors, sarParameters.tstart, tstop);
else
    XPoly = sarParameters.satelliteTrajectory.XPoly;
end
XSat = polyvalSV(XPoly, azimuthTime); 

% Line of sight vector
LOS = PointXYZ - XSat;

% Complementary angle of off-local vertical (incidence) angle
cosIncAngle = dot(PointXYZ,LOS);
sinIncAngle = norm(cross(PointXYZ,LOS));
incAngle = pi - atan2(sinIncAngle,cosIncAngle);           
theta = pi/2 - incAngle;

% Compute coordinates of two close by lines on ground
[latP, lonP] = CalcSARLatLon(sample, line - 1, height, sarParameters);
[xP, yP, zP] = ell2xyz(deg2rad(latP), deg2rad(lonP), height, elliSemiMajor, elliEccSquared);
Point1XYZ = [xP; yP; zP];

[latP, lonP] = CalcSARLatLon(sample, line + 1, height, sarParameters);
[xP, yP, zP] = ell2xyz(deg2rad(latP), deg2rad(lonP), height, elliSemiMajor, elliEccSquared);
Point2XYZ = [xP; yP; zP];

% Compute heading
N = [0; 0; 1]; % Z unit vector
% normal to great circle passing through P1 and P2 
c1 = cross(Point2XYZ, Point1XYZ);
% normal to great circle passing through P1 and N 
c2 = cross(Point1XYZ, N);
sinHeading = norm(cross(c1,c2));
cosHeading = dot(c1,c2);

heading = atan2(sinHeading,cosHeading) - pi;
phi = -(pi + heading);

% % Wrap in [-pi,pi]
% if phi >= 0
%     phi = phi - 2*pi*floor((phi + pi) / (2*pi))
% else
%     phi = phi - 2*pi*ceil ((phi - pi) / (2*pi));
% end

