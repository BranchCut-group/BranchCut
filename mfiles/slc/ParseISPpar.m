function [S] = ParseISPpar(slcPar)
%
% Creator John Merryman - INGV
% Date: 02 May 2013
%
% Usage: S = ParseISPpar(SLC_PAR)
%
% Parse GAMMA ISP (Interferometric SAR processor) Single Look Complex
% ASCII parameter file.
%
% SLC_PAR   : (input)  GAMMA ISP SLC ASCII parameter file path.
% S         : (output) Struct containing relevant SLC parameters.
%
% Example:
%
% S = ParseISPpar('../data/19960410.rslc.par');
%

fprintf('\nReading file: %s\n', slcPar);

fid = fopen(slcPar);

A = textscan(fid, '%s%s', 'Delimiter', ':', 'Headerlines', 1);

fclose(fid);

nli = size(A{1,1}, 1);

S = struct;
S.fn = slcPar;

fprintf('\nStructure dump:\n\n');
for i = 1:nli

    if strcmp(A{1,1}{i,1}, 'sensor')
        token = strtok(A{1,2}{i,1});
        S.sensor = token;
        fprintf('    sensor: %s \n', S.sensor);
        % Assign antenna length
        if strcmp(S.sensor, 'CSKS1') || strcmp(S.sensor, 'CSKS2') || ...
           strcmp(S.sensor, 'CSKS3') || strcmp(S.sensor, 'CSKS4') || ...
           strcmp(S.sensor, 'COSMO-SkyMed')
            S.L = 5.7;
        elseif strcmp(S.sensor, 'S1A')
            S.L = 12.3;
        else
            S.L = 10;
        end
        fprintf('    L: %f m\n', S.L);
    end
    
    if strcmp(A{1,1}{i,1}, 'date')
        token = strtok(A{1,2}{i,1}, '\n');
        yyyymmdd = str2num(token);
        dateStr = upper(datestr([yyyymmdd(1) yyyymmdd(2) yyyymmdd(3) 0 0 0]));
        dateSec = date2sec(dateStr);
        fprintf('    dateStr = %s\n', dateStr);
    end     
    
    if strcmp(A{1,1}{i,1}, 'start_time')
        token = strtok(A{1,2}{i,1}, 's');
        S.tstart = dateSec + str2num(token);
        fprintf('    tstart: %f s\n', S.tstart);
    end
    
    if strcmp(A{1,1}{i,1}, 'center_time')
        token = strtok(A{1,2}{i,1}, 's');
        S.tcenter = dateSec + str2num(token);
        fprintf('    tcenter: %f s\n', S.tcenter);
    end
    
    if strcmp(A{1,1}{i,1}, 'end_time')
        token = strtok(A{1,2}{i,1}, 's');
        S.tstop = dateSec + str2num(token);
        fprintf('    tstop: %f s\n', S.tstop);
    end 
    if strcmp(A{1,1}{i,1}, 'azimuth_line_time')
        token = strtok(A{1,2}{i,1}, 's');
        S.PRT = str2num(token);
        fprintf('    PRT: %e s\n', S.PRT);
    end 

    if strcmp(A{1,1}{i,1}, 'range_samples')
        token = strtok(A{1,2}{i,1});
        S.nsa = str2num(token);
        fprintf('    nsa: %d \n', S.nsa);
    end
    
    if strcmp(A{1,1}{i,1}, 'azimuth_lines')
        token = strtok(A{1,2}{i,1});
        S.nli = str2num(token);
        fprintf('    nli: %d \n', S.nli);
    end 
    
    if strcmp(A{1,1}{i,1}, 'radar_frequency')
        token = strtok(A{1,2}{i,1}, 'Hz');
        S.f0 = str2num(token);
        SOL = 299792458;
        S.SOL = SOL;
        S.lambda = SOL / S.f0;
        fprintf('    f0: %e Hz\n', S.f0);
        fprintf('    lambda: %f m\n', S.lambda);
    end 

    if strcmp(A{1,1}{i,1}, 'center_latitude')
        token = strtok(A{1,2}{i,1}, 'degrees');
        S.lat0 = str2num(token);
        fprintf('    lat0: %f deg\n', S.lat0);
    end     

    if strcmp(A{1,1}{i,1}, 'center_longitude')
        token = strtok(A{1,2}{i,1}, 'degrees');
        S.lon0 = str2num(token);
        fprintf('    lon0: %f deg\n', S.lon0);
    end     

    if strcmp(A{1,1}{i,1}, 'heading')
        token = strtok(A{1,2}{i,1}, 'degrees');
        S.heading = str2num(token);
        fprintf('    heading: %f deg\n', S.heading);
    end 
    
    if strcmp(A{1,1}{i,1}, 'prf')
        token = strtok(A{1,2}{i,1}, 'Hz');
        S.PRF = str2num(token);
        fprintf('    PRF: %f Hz\n', S.PRF);
    end 

    if strcmp(A{1,1}{i,1}, 'near_range_slc')
        token = strtok(A{1,2}{i,1}, 'm');
        S.Rnear = str2num(token);
        fprintf('    Rnear: %f m\n', S.Rnear);
    end 

    if strcmp(A{1,1}{i,1}, 'center_range_slc')
        token = strtok(A{1,2}{i,1}, 'm');
        S.Rcen = str2num(token);
        fprintf('    Rcen: %f m\n', S.Rcen);
    end 

    if strcmp(A{1,1}{i,1}, 'far_range_slc')
        token = strtok(A{1,2}{i,1}, 'm');
        S.Rfar = str2num(token);
        fprintf('    Rfar: %f\n', S.Rfar);
    end 

    if strcmp(A{1,1}{i,1}, 'incidence_angle')
        token = strtok(A{1,2}{i,1}, 'degrees');
        S.inc_angle = str2num(token);
        fprintf('    incidence_angle: %f deg\n', S.inc_angle);
    end 
    
    if strcmp(A{1,1}{i,1}, 'azimuth_proc_bandwidth')
        token = strtok(A{1,2}{i,1}, 'Hz');
        S.Bp = str2num(token);
        fprintf('    Bp: %f Hz\n', S.Bp);
    end 
    
    if strcmp(A{1,1}{i,1}, 'adc_sampling_rate')
        token = strtok(A{1,2}{i,1}, 'Hz');
        S.fs = str2num(token);
        fprintf('    fs: %f Hz\n', S.fs);
    end     
    
    if strcmp(A{1,1}{i,1}, 'chirp_bandwidth')
        token = strtok(A{1,2}{i,1}, 'Hz');
        S.Br = str2num(token);
        fprintf('    Br: %f Hz\n', S.Br);
    end    
    
    if strcmp(A{1,1}{i,1}, 'number_of_state_vectors')
        token = strtok(A{1,2}{i,1});
        S.Nsv = str2num(token);
        fprintf('    Nsv: %d \n', S.Nsv);
    end 

    if strcmp(A{1,1}{i,1}, 'time_of_first_state_vector')
        token = strtok(A{1,2}{i,1}, 's');
        S.t0sv = dateSec + str2num(token);
        fprintf('    t0sv: %f s\n', S.t0sv);
    end 
    
    if strcmp(A{1,1}{i,1}, 'state_vector_interval')
        token = strtok(A{1,2}{i,1}, 's');
        S.dtsv = str2num(token);
        fprintf('    dtsv: %f s\n', S.dtsv);
    end 

    if strcmp(A{1,1}{i,1}, 'doppler_polynomial')
        token = strtok(A{1,2}{i,1}, 'Hz');
        S.dopplerCentroidPoly.rangeCoeffs = str2num(token);
        S.dopplerCentroidPoly.t0 = 0;
        fprintf('    fdcp: %f Hz; %f Hz/m; %f Hz/m^2;\n', ...
                 S.dopplerCentroidPoly.rangeCoeffs(1), S.dopplerCentroidPoly.rangeCoeffs(2), ...
                 S.dopplerCentroidPoly.rangeCoeffs(3));
    end     

    if strcmp(A{1,1}{i,1}, 'range_pixel_spacing')
        token = strtok(A{1,2}{i,1}, 'm');
        S.dr = str2num(token);
        fprintf('    dr: %f m\n', S.dr);
    end

    if strcmp(A{1,1}{i,1}, 'azimuth_pixel_spacing')
        token = strtok(A{1,2}{i,1}, 'm');
        S.dx = str2num(token);
        fprintf('    dx: %f m\n', S.dx);
    end

    if strcmp(A{1,1}{i,1}, 'sar_to_earth_center')
        token = strtok(A{1,2}{i,1}, 'm');
        S.H = str2num(token);
        fprintf('    H: %f m\n', S.H);
    end    

    if strcmp(A{1,1}{i,1}, 'earth_radius_below_sensor')
        token = strtok(A{1,2}{i,1}, 'm');
        S.Re = str2num(token);
        fprintf('    Re: %f m\n', S.Re);
    end
    
    if strcmp(A{1,1}{i,1}, 'azimuth_angle')
        token = strtok(A{1,2}{i,1}, 'degrees');
        if str2num(token) == 90;
            S.lookDirection = 'RIGHT';
        else
            S.lookDirection = 'LEFT';            
        end
        fprintf('    look direction: %s \n', S.lookDirection);
    end     
    
end
fprintf('\n');

for i = 1:nli
     for j = 1:S.Nsv
         S.stateVectors.tsv(j,1) = S.t0sv + (j - 1)*S.dtsv;
         param = sprintf('state_vector_position_%d',j);    
         if strcmp(A{1,1}{i,1}, param)
             token = strtok(A{1,2}{i,1}, 'm');
             S.stateVectors.X(j,:) = str2num(token);
         end 
         param = sprintf('state_vector_velocity_%d',j);    
         if strcmp(A{1,1}{i,1}, param)
             token = strtok(A{1,2}{i,1}, 'm');
             S.stateVectors.V(j,:) = str2num(token);
         end 
     end
end
