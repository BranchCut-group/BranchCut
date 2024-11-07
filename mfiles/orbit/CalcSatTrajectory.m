function [XPoly, VPoly] = CalcSatTrajectory(stateVectors, tstart, tstop, ...
    polyDegree, maxTimeMargin, doScale)
%
% Creator John Peter Merryman Boncori - sarmap 
% Date: 24 Jan 2017
%
% usage: [XPoly, VPoly] = CalcSatTrajectory(stateVectors, tstart, tstop, ...
%     polyDegree, maxTimeMargin)
%
% Fit polynomials to compute satellite position and velocity as a
% function of azimuth time.
%
% stateVectors   : (input)  Satellite state vector struct
% tstart         : (input)  Data start time 
% tstop          : (input)  Data stop time 
% polyDegree     : (input)  degree of polynomial fit (default = 4)
% maxTimeMargin  : (input)  Time interval in seconds. Only state vectors
%                           within maxTimeMargin from tstart and tstop are used
% doScale        : (input)   centre and scale variable before fitting
% XPoly          : (output) Position polynomial coefficients
% VPoly          : (output) Velocity poynomial coefficients
%
% To evaluate position and velocity at given time tazi (in seconds of day):
%   Xsat = polyvalSV(XPoly, tazi); 
%   Vsat = polyvalSV(VPoly, tazi);
%
  if nargin < 3
      error('Insufficient number of input parameters, exiting.');
  end

  if nargin < 4
      polyDegree = 4;
  end

  if nargin < 5
      maxTimeMargin = 120;
  end    

  if nargin < 6
     doScale = 1;
  end  
  
  tmid = (tstart + tstop) / 2;
  dt = tstop - tstart;
  
  svStep = stateVectors.tsv(2) - stateVectors.tsv(1);
  maxTimeMargin = max(maxTimeMargin, (polyDegree + 2) / 2 * svStep);
  
  nearbyStateVectorIdx = find(abs(stateVectors.tsv - tmid) <= dt/2 + maxTimeMargin);
  numOfNearbyStateVectors = length(nearbyStateVectorIdx);
  if numOfNearbyStateVectors < polyDegree + 1
     error('Insufficient number of state vectors (= %d), exiting.', numOfNearbyStateVectors);
  end 

  nearbyStateVectors.tsv = stateVectors.tsv(nearbyStateVectorIdx, 1);
  nearbyStateVectors.X   = stateVectors.X(nearbyStateVectorIdx, :);
  nearbyStateVectors.V   = stateVectors.V(nearbyStateVectorIdx, :);
  
%   for i = 1 : numOfNearbyStateVectors
%       fprintf('\nState vector %d\n', nearbyStateVectorIdx(i))
%       fprintf('t  = %.20f\n', nearbyStateVectors.tsv(i) - stateVectors.tsv(1));
%       fprintf('x  = %.20f\n', nearbyStateVectors.X(i,1));
%       fprintf('y  = %.20f\n', nearbyStateVectors.X(i,2));
%       fprintf('z  = %.20f\n', nearbyStateVectors.X(i,3));
%       fprintf('vx = %.20f\n', nearbyStateVectors.V(i,1));
%       fprintf('vy = %.20f\n', nearbyStateVectors.V(i,2));
%       fprintf('vz = %.20f\n', nearbyStateVectors.V(i,3));
%   end
  
  % Fit satellite position and velocity as a function of azimuth time
  XPoly = polyfitSV(nearbyStateVectors.X, nearbyStateVectors.tsv, polyDegree, doScale);
  VPoly = polyfitSV(nearbyStateVectors.V, nearbyStateVectors.tsv, polyDegree, doScale);
