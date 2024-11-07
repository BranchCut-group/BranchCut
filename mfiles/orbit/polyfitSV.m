function  P = polyfitSV(X, t, polyDegree, doScale)
%
% Creator John Peter Merryman Boncori - sarmap 
% Date: 14 Jul 2016
%
% usage: P = polyfitSV(X, t, polyDegree);
%
% Polynomial fit of vector components.
%
% X          : (input)   state vectors: Nsv x 3 matrix. Nsv is number of state
%                        vectors. Each column contains an x,y,z component.
% t          : (input)   state vector times
% polyDegree : (input)   degree of polynomial fit
% doScale    : (input)   centre and scale variable before fitting
% P          : (output)  polyfit struct 
%
% Example:
%     XPoly = polyfitSV(stateVectors.X, stateVectors.tsv, polyDegree);
%

if nargin < 3
    polyDegree = 4;
end

if nargin < 4
    doScale = 1;
end

% Polynomial fit of state vector info with t as independent variable
if doScale == 1
    [P.px, P.Sx, P.mux] = polyfit(t, X(:,1), polyDegree);
    [P.py, P.Sy, P.muy] = polyfit(t, X(:,2), polyDegree);
    [P.pz, P.Sz, P.muz] = polyfit(t, X(:,3), polyDegree);
else
    [P.px, P.Sx] = polyfit(t, X(:,1), polyDegree);
    [P.py, P.Sy] = polyfit(t, X(:,2), polyDegree);
    [P.pz, P.Sz] = polyfit(t, X(:,3), polyDegree);
    P.mux = [0 1];
    P.muy = [0 1];
    P.muz = [0 1];    
end    