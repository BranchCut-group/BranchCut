function [S] = ParseDEMpar(demPar)
%
% Creator John Merryman - INGV
% Date: 12 Jan 2015
%
% Usage: [S] = ParseDEMpar (DEM_PAR)
%
% Parse GAMMA DIFF&GEO DEM ASCII parameter file.
%
% DEM_PAR   : (input)  GAMMA DIFF&GEO DEM ASCII parameter file
% S         : (output) Struct containing relevant parameters.
%

fprintf('\nReading file: %s\n', demPar);

fid = fopen(demPar);

A = textscan(fid, '%s%s', 'Delimiter', ':', 'Headerlines', 1);

fclose(fid);

nli = size(A{1,1}, 1);

S = struct;
S.fn = demPar;

fprintf('\nStructure dump:\n\n');
for i = 1:nli    
    
    if strcmp(A{1,1}{i,1}, 'nlines')
        token = strtok(A{1,2}{i,1});
        S.nrows = str2num(token);
        fprintf('    nrows: %d \n', S.nrows);
    end
    
    if strcmp(A{1,1}{i,1}, 'width')
        token = strtok(A{1,2}{i,1});
        S.ncols = str2num(token);
        fprintf('    ncols: %d \n', S.ncols);
    end

    if strcmp(A{1,1}{i,1}, 'corner_lat')
        token = strtok(A{1,2}{i,1});
        S.corner_lat = str2num(token);
        fprintf('    corner_lat: %f decimal degrees\n', S.corner_lat);
    end    

    if strcmp(A{1,1}{i,1}, 'corner_lon')
        token = strtok(A{1,2}{i,1});
        S.corner_lon = str2num(token);
        fprintf('    corner_lon: %f decimal degrees\n', S.corner_lon);
    end
    
    if strcmp(A{1,1}{i,1}, 'post_lat')
        token = strtok(A{1,2}{i,1});
        S.post_lat = str2num(token);
        fprintf('    post_lat: %e decimal degrees\n', S.post_lat);
    end    
    
    if strcmp(A{1,1}{i,1}, 'post_lon')
        token = strtok(A{1,2}{i,1});
        S.post_lon = str2num(token);
        fprintf('    post_lon: %e decimal degrees\n', S.post_lon);
    end    
    
end