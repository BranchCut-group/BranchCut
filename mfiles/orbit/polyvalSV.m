function  X = polyvalSV(P, t0, doScale)
%
% Creator John Peter Merryman Boncori - sarmap 
% Date: 14 Jul 2016
%
% usage: X = polyvalSV(P, t0);
%
% Evaluate vector component polynomial coefficients at a given time.
%
% P       : (input)  polyfit struct
% t0      : (input)  independent variable value for which polynomial fit must be evaluated
% doScale : (input)  centre and scale variable before fitting
% X       : (output) vector components corresponding to t0
%
% Example:
%     X = polyvalSV(XPoly, t0); 
%
%

if nargin < 3
    doScale = 1;
end

if doScale == 1
  x = polyval(P.px, t0, P.Sx, P.mux);
  y = polyval(P.py, t0, P.Sy, P.muy);
  z = polyval(P.pz, t0, P.Sz, P.muz);
else
  x = polyval(P.px, t0);
  y = polyval(P.py, t0);
  z = polyval(P.pz, t0);    
end

X = [x; y; z];
