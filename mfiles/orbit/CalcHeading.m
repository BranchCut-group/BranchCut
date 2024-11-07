function heading = CalcHeading(sarParameters, line)
%
% Creator John Peter Merryman Boncori - sarmap 
% Date: 14 Feb 2017
%
% Usage: 
%
% heading = CalcHeading(sarParameters, line);
%
% Compute heading angle for a azimuth line.
% The heading angle is counted from the north direction and should be 
% about -15 deg for ascending orbits and -165 deg for descending ones.
% 
% sarParameters     : (input)  Raw (or SLC) data parameter struct
% line              : (input)  Raw or SLC azimuth line number
% heading           : (output) Heading angle (deg)
%
% Example:
%     S = ParseRawSentinelSMLpar('sentinel_queensland_IW2_5.sml');
%     heading = CalcHeading(S, S.nli/2);
%

if nargin < 2
  error('Insufficient number of input parameters, exiting.');
end

% Compute coordinates of two close by points along the orbit
azimuthTime = sarParameters.tstart + (line - 1) / sarParameters.PRF;
dt = 0.01;
Point2XYZ = polyvalSV(sarParameters.satelliteTrajectory.XPoly, azimuthTime + dt); 
Point1XYZ = polyvalSV(sarParameters.satelliteTrajectory.XPoly, azimuthTime - dt); 

% Compute heading
N = [0; 0; 1]; % Z unit vector
% normal to great circle passing through P1 and P2 
c1 = cross(Point2XYZ, Point1XYZ);
% normal to great circle passing through P1 and N 
c2 = cross(Point1XYZ, N);
sinHeading = norm(cross(c1,c2));
cosHeading = dot(c1,c2);

heading = rad2deg(atan2(sinHeading,cosHeading) - pi);
