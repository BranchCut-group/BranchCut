function geocodedImage = geocodeBack(radarImage, geo2SarLUT)
%
% Creator John Peter Merryman Boncori - DTU 
% Date: 28 May 2018
%
% Usage:
%     geocodedImage = geocodeBack(radarImage, geo2SarLUT)
%
% Geocode image in radar geometry based on lookup table in geocoded geometry
%
% radarImage     : (input)  Image in radar geometry
% LUT            : (input)  Complex lookup table in geocoded geometry.
% geocodedImage  : (output) Image in geocoded geometry.
%

geocodedImage = interp2(radarImage, real(geo2SarLUT), imag(geo2SarLUT));
