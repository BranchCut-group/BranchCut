function wrappedPhase = wrap(phase)
%
% Wrap phase matrix in [-pi, pi).
%
% Creator John Merryman - DTU
% Date: 25 Oct. 2017
%
% Usage:
%  w = wrap(pha);
%
% phase         :  (input)  Phase matrix (rad)
% wrappedPhase  : (output)  Wrapped phase matrix (rad).
%

wrappedPhase = mod(phase + pi, 2*pi) - pi;
