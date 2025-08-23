close all; clear; clc;

eph = read_GPSbroadcast("data_gps" + filesep + "brdc2930.11n");
eph = sortrows(eph, 1);

prn_list = [2, 4, 5, 9, 10, 12, 17, 23, 25];

for ii = 1:numel(prn_list)
    prn = prn_list(ii);
    slice = eph(eph(:, 1) == prn, :);
    gps_week = unique(slice(:, 19));
    
    temp = transpose(417088:1:417986);
    gps_week_vec = gps_week * ones(size(temp));
    time_vec = [gps_week_vec, temp];
    t_steps = size(time_vec, 1);
    
    [health, x, v, relcorr, satClkCorr] = broadcast2xv(eph, time_vec, prn, 'eph');
    
    save_data = nan(size(x, 1), 12);
    save_data(:, 1) = prn * ones(t_steps, 1);
    save_data(:, 2) = gps_week * ones(t_steps, 1);
    save_data(:, 3) = temp;
    save_data(:, 4:6) = x;
    save_data(:, 7:9) = v;
    save_data(:, 10) = relcorr;
    save_data(:, 11) = satClkCorr;
    save_data(:, 12) = health;

    save_data = array2table(save_data,...
        'VariableNames', {'PRN', 'GPS_WEEK', 'TOC',...
        'SAT_POS_X_M', 'SAT_POS_Y_M', 'SAT_POS_Z_M',...
        'SAT_VEL_X_MPS', 'SAT_VEL_Y_MPS', 'SAT_VEL_Z_MPS',...
        'REL_CLK_CORR', 'SAT_CLK_CORR', 'HEALTH'});

    writetable(save_data, "data_gps" + filesep +...
        "gps_sat_state_data_prn_" + string(prn) + "_.csv");
    
end


% for jj = 1:numel(prn_list)
%     prn = prn_list(jj);
%     slice = eph(eph(:, 1) == prn, :);
%     gps_week = unique(slice(:, 19));
%     gps_time = unique(slice(:, 20));
%     for ii = 1:(numel(gps_time)-1)
%         temp = transpose(gps_time(ii):1:gps_time(ii+1));
%         gps_week_vec = gps_week * ones(size(temp));
%         time_vec = [gps_week_vec, temp];
%         t_steps = size(time_vec, 1);
% 
%         [health, x, v, relcorr, satClkCorr] = broadcast2xv(eph, time_vec, prn, 'eph');
% 
%         save_data = nan(size(x, 1), 12);
%         save_data(:, 1) = prn * ones(t_steps, 1);
%         save_data(:, 2) = gps_week * ones(t_steps, 1);
%         save_data(:, 3) = temp;
%         save_data(:, 4:6) = x;
%         save_data(:, 7:9) = v;
%         save_data(:, 10) = relcorr;
%         save_data(:, 11) = satClkCorr;
%         save_data(:, 12) = health;
% 
%         save_data = array2table(save_data,...
%             'VariableNames', {'PRN', 'GPS_WEEK', 'TOC',...
%             'SAT_POS_X_M', 'SAT_POS_Y_M', 'SAT_POS_Z_M',...
%             'SAT_VEL_X_MPS', 'SAT_VEL_Y_MPS', 'SAT_VEL_Z_MPS',...
%             'REL_CLK_CORR', 'SAT_CLK_CORR', 'HEALTH'});
% 
%         writetable(save_data, "data_gps" + filesep +...
%             "gps_sat_state_data_prn" + string(prn) +...
%             "_gps_week_" +  string(gps_week) +...
%             "_start_" + string(gps_time(ii)) +...
%             "_end_" + string(gps_time(ii+1)) + "_.csv");
% 
%     end
% end

%% Save to csv for porting over to python
eph = array2table(eph,...
    'VariableNames', {'PRN', 'M0', 'DELTA_N', 'ECC', 'SQRT_A', 'LOA',...
    'INC', 'PERIGEE', 'RA_RATE', 'INC_RATE', 'Cuc', 'Cus', 'Crc', 'Crs',...
    'Cic', 'Cis', 'Toe', 'IODE', 'GPS_WEEK', 'TOC', 'Af0', 'Af1', 'Af2',...
    'BLANK', 'HEALTH'});
writetable(eph, "data_gps" + filesep + "gps_broadcast_data.csv");

%% Helper Functions
function [gps_ephem,ionoparams] = read_GPSbroadcast(navfilename)
%==========================================================================
%==========================================================================
% [gps_ephem,ionoparams] = read_GPSbroadcast(navfilename)
%
% Reads an IGS GPS Broadcast Ephemeris for all GPS satellites and 
%  constructs a matrix of all ephemeris values. Note that each row of the 
%  gps_ephem is a new entry in the broadcast file which will most likely 
%  contain multiple entries for the same PRN.
%
%
% Author: Ben K. Bradley
% Date: 07/19/2009
%
%
% INPUT:             Description                                      Units
%
%  navfilename   - name of IGS broadcast ephemeris file to read in   string
%
%
% OUTPUT:       
%    
%  gps_ephem     - matrix of gps satellite orbit parameters          (nx25)
%  
%                  col1: prn, PRN number of satellite
%                  col2: M0, mean anomaly at reference time, rad
%                  col3: delta_n, mean motion difference from computed value, rad/s
%                  col4: ecc, eccentricity of orbit
%                  col5: sqrt_a, square root of semi-major axis, m^0.5
%                  col6: Loa, longitude of ascending node of orbit plane at weekly epoch, rad
%                  col7: incl, inclination angle at reference time, rad
%                  col8: perigee, argument of perigee, rad
%                  col9: ra_rate, rate of change of right ascension, rad/s
%                 col10: i_rate, rate of change of inclination angle, rad/s
%                 col11: Cuc, amplitude of the cosine harmonic correction term to the argument of latitude
%                 col12: Cus, amplitude of the sine harmonic correction term to the argument of latitude
%                 col13: Crc, amplitude of the cosine harmonic correction term to the orbit radius
%                 col14: Crs, amplitude of the sine harmonic correction term to the orbit radius
%                 col15: Cic, amplitude of the cosine harmonic correction term to the angle of inclination
%                 col16: Cis, amplitude of the cosine harmonic correction term to the angle of inclination
%                 col17: Toe, reference time ephemeris (seconds into GPS week)
%                 col18: IODE, issue of data (ephemeris) 
%                 col19: GPS_week, GPS Week Number (to go with Toe) wrt 06jan80 (no rollover)
%                 col20: Toc, time of clock
%                 col21: Af0, satellite clock bias (sec)
%                 col22: Af1, satellite clock drift (sec/sec)
%                 col23: Af2, satellite clock drift rate (sec/sec/sec)
%                 col24: blank (zero)
%                 col25: health, satellite health (0=good and usable)
%
%
%  ionoparams    - parameters for the Klobuchar  [A0 A1 A2 A3 B0 B1 B2 B3]
%                    ionospheric model                        
%          
%
% Coupling:
%
%  none
%
% References:
% 
%  [1] Interface Control Document: IS-GPS-200D
%        < http://www.navcen.uscg.gov/gps/geninfo/IS-GPS-200D.pdf >
%
%  [2] RINEX GPS Format, Version 2, (Table A4)
%        < http://www.ngs.noaa.gov/CORS/instructions2/ >
%
%==========================================================================
%==========================================================================

% Open desired IGS Ephemeris File
%==========================================================================

if (exist(navfilename,'file') == 2)

    fid = fopen(navfilename,'r');
else
    error(sprintf('Unable to find broadcast file: %s',navfilename), 'ERROR!');
end


% Step through the header of the file and pull out iono parameters
%==========================================================================
headerend   = [];
headeralpha = [];  ALPHA = [];
headerbeta  = [];  BETA  = [];

while (isempty(headerend) == 1)
   tline     = fgetl(fid); 
   headerend = findstr(tline,'END OF HEADER');
   
   headeralpha = findstr(tline,'ION ALPHA');
   if (isempty(headeralpha) == 0)
       
      [A0, remain] = strtok(tline);
      [A1, remain] = strtok(remain);
      [A2, remain] = strtok(remain);
      [A3]         = strtok(remain);
      
      ALPHA = [str2num(A0) str2num(A1) str2num(A2) str2num(A3)];
      
   end
   
   headerbeta = findstr(tline,'ION BETA');
   if (isempty(headerbeta) == 0)
       
      [B0, remain] = strtok(tline);
      [B1, remain] = strtok(remain);
      [B2, remain] = strtok(remain);
      [B3]         = strtok(remain);
      
      BETA = [str2num(B0) str2num(B1) str2num(B2) str2num(B3)];
      
   end
   
end

ionoparams = [ALPHA BETA];

j = 1;

% Enter main loop to read the rest of the ephemeris file
%==========================================================================
%==========================================================================
while 1
    
    % Load next line in ephemeris file
    tline = fgetl(fid);
    
    % If the next line is not a character then the end of the file has been
    %   reached and the while loop is exited
    if ~ischar(tline), break, end
   
        %-----------------------------------------------------------------
        % Read in variables of the FIRST line of this satellite's ephemeris
        %-----------------------------------------------------------------
        prn = str2num(tline(1:2));   % PRN number of satellite
       
        Af0 = str2num(tline(23:41)); % clock bias  (s)
        Af1 = str2num(tline(42:60)); % clock drift (s/s)
        Af2 = str2num(tline(61:79)); % clock drift rate (s/s/s)      
        
        %-----------------------------------------------------------------
        % SECOND LINE
        %-----------------------------------------------------------------
        tline = fgetl(fid); % Read in second line of satellite ephemeris
        
        IODE = str2num(tline(4:22));    % Issue of Data (Ephemeris)  
        
        Crs  = str2num(tline(23:41));   % Amplitude of the Sine Harmonic Correction 
                                        %  Term to the orbit radius
      
        delta_n= str2num(tline(42:60)); % Mean Motion Difference from Computed Value, rad/s
        
        M0   = str2num(tline(61:79));   % Mean Anomaly at Reference Time, rad
        
        %-----------------------------------------------------------------
        % THIRD LINE
        %-----------------------------------------------------------------
        tline = fgetl(fid); % Read in third line of satellite ephemeris
        
        Cuc = str2num(tline(4:22));     % Amplitude of the Cosine Harmonic Correction
                                        %  Term to the Argument of Latitude
       
        ecc = str2num(tline(23:41));    % Eccentricity
        
        Cus = str2num(tline(42:60));    % Amplitude of the Sine Harmonic Correction
                                        %  Term to the Argument of Latitude
       
        sqrt_a = str2num(tline(61:79)); % Square root of the semi-major Axis, m^0.5
        
        %-----------------------------------------------------------------
        % FOURTH LINE
        %-----------------------------------------------------------------
        tline = fgetl(fid); % Read in fourth line of satellite ephemeris
        
        Toe = str2num(tline(4:22));     % Reference Time Ephemeris (sec into GPS week)
        
        Cic = str2num(tline(23:41));    % Amplitude of the Cosine Harmonic Correction
                                        %  Term to the Angle of Inclination
        
        Loa = str2num(tline(42:60));    % Longitude of Ascending Node of Orbit Plane
                                        %  at Weekly Epoch, rad
      
        Cis = str2num(tline(61:79));    % Amplitude of the Sine Harmonic Correction
                                        %  Term to the Angle of Inclination
         
        %-----------------------------------------------------------------                    
        % FIFTH LINE
        %-----------------------------------------------------------------
        tline = fgetl(fid); % Read in fifth line of satellite ephemeris
        
        incl = str2num(tline(4:22));    % Inclination Angle at Reference Time, rad
        
        Crc  = str2num(tline(23:41));   % Amplitude of the Cosine Harmonic Correction
                                        %  Term to the Orbit Radius
     
        perigee = str2num(tline(42:60));% Argument of Perigee, rad
        
        ra_rate = str2num(tline(61:79));% Rate of Change of Right Ascension, rad/s
        
        %-----------------------------------------------------------------
        % SIXTH LINE
        %-----------------------------------------------------------------
        tline = fgetl(fid); % Read in sixth line of satellite ephemeris
        
        i_rate = str2num(tline(4:22));   % Rate of change of inclination angle, rad/s
        
        %str = tline(23:41);             % codes on L2 channel (unecessary)
       
        GPS_week = str2num(tline(42:60));% GPS Week Number (to go with Toe)
        if GPS_week < 1024
            GPS_week = GPS_week+1024; % Penny added to avoid modulo 1024
        end
        %str   = tline(61:79);           % L2 flag
        
        %-----------------------------------------------------------------
        % SEVENTH LINE 
        %-----------------------------------------------------------------
        tline = fgetl(fid); % Includes: SV accuracy, SV health, TGD, IODC
        
        health = str2num(tline(23:41)); % Satellite health (0.00 = usable)
        
        %-----------------------------------------------------------------
        % EIGHTH LINE 
        %-----------------------------------------------------------------
        tline = fgetl(fid); % Read in eighth line of satellite ephemeris
        
        Toc = Toe; % Time of clock   
        
        gps_ephem(j,:) = [prn M0 delta_n ecc sqrt_a Loa incl perigee ra_rate i_rate Cuc Cus Crc Crs Cic Cis Toe IODE GPS_week Toc Af0 Af1 Af2 0 health];
        
        j = j + 1;   
        
end


fclose(fid);
end

function [health,x,v,relcorr,satClkCorr] = broadcast2xv(ephem_all,t_input,prn, input)

%==========================================================================
%==========================================================================
% [health,x,v,relcorr,satClkCorr] = broadcast2xv(ephem_all,t_input,prn)
%
% Calculates the position and velocity of a GPS satellite from an ephemeris 
%  matrix (see read_GPSbroadcast.m).  The relativity correction and 
%  satellite clock correction are also computed.  The input ephem_all can 
%  be generated by the read_GPSbroadcast.m function.
%
% When computing an expected pseudorange from receiver to GPS satellite,
%  the following sign convention should be used when using the relativity
%  and satellite clock corrections:
%
%  pseudorange = RANGE - SATclkCORR + REL_CORR + ...
%
%
%
% Author: Ben K. Bradley
% Date: 07/19/2009
% Mod: Sara Hrbek 
% Date: 9/24/2016
% INPUT:               Description                                  Units
%
%  ephem_all    - matrix of gps satellite orbit parameters           (nx25)
%  
%                  col1: prn, PRN number of satellite
%                  col2: M0, mean anomaly at reference time, rad
%                  col3: delta_n, mean motion difference from computed value, rad/s
%                  col4: ecc, eccentricity of orbit
%                  col5: sqrt_a, square root of semi-major axis, m^0.5
%                  col6: Loa, longitude of ascending node of orbit plane at weekly epoch, rad
%                  col7: incl, inclination angle at reference time, rad
%                  col8: perigee, argument of perigee, rad
%                  col9: ra_rate, rate of change of right ascension, rad/s
%                 col10: i_rate, rate of change of inclination angle, rad/s
%                 col11: Cuc, amplitude of the cosine harmonic correction term to the argument of latitude
%                 col12: Cus, amplitude of the sine harmonic correction term to the argument of latitude
%                 col13: Crc, amplitude of the cosine harmonic correction term to the orbit radius
%                 col14: Crs, amplitude of the sine harmonic correction term to the orbit radius
%                 col15: Cic, amplitude of the cosine harmonic correction term to the angle of inclination
%                 col16: Cis, amplitude of the cosine harmonic correction term to the angle of inclination
%                 col17: Toe, reference time ephemeris (seconds into GPS week)
%                 col18: IODE, issue of data (ephemeris) 
%                 col19: GPS_week, GPS Week Number (to go with Toe)
%                 col20: Toc, time of clock
%                 col21: Af0, satellite clock bias (sec)
%                 col22: Af1, satellite clock drift (sec/sec)
%                 col23: Af2, satellite clock drift rate (sec/sec/sec)
%                 col24: blank (zero)
%                 col25: health, satellite health (0=good and usable)
%
%  t_input      - GPS times to calculate values at                 [WN TOW] (nx2)
%  prn          - PRN to compute values for (one satellite only)                       
%  input        - 'alm' for almanac or 'eph' for ephemeris
%
%
% OUTPUT:       
%    
%  health       - health of satellite (0=good)                              (nx1)
%  x            - position of satellite (ECEF)                  [x y z]   m (nx3)
%  v            - velocity of satellite (ECEF, frame & mag)  [dx dy dz] m/s (nx3)
%  relcorr      - relativity correction                                   m (nx1)
%  satClkCorr   - GPS satellite clock correction (w/o relcorr)            m (nx1)
%                                     
%
%
% Coupling:
%
%   mean2eccentric.m
%
% References:
% 
%   [1] Interface Control Document: IS-GPS-200D
%         < http://www.navcen.uscg.gov/gps/geninfo/IS-GPS-200D.pdf >
%
%   [2] Zhang, J., et.all. "GPS Satellite Velocity and Acceleration
%         Determination using the Broadcast Ephemeris". The Journal of
%         Navigation. (2006), 59, 293-305.
%            < http://journals.cambridge.org/action/displayAbstract;jsess ...
%                ionid=C6B8C16A69DD7C910989C661BAB15E07.tomcat1?fromPage=online&aid=425362 >
%
%   [3] skyplot.cpp by the National Geodetic Survey
%          < http://www.ngs.noaa.gov/gps-toolbox/skyplot/skyplot.cpp >
%
%
% Last Updated:
%
%  2015/01/22  B.K. Bradley - the capability to look for updated ephem
%                              entries that occur at odd times within each
%                              2hr window has been commented out in this 
%                              function and added to read_GPSbroadcast.m
%                              instead. This moves the computational
%                              overhead to the reading which only occurs
%                              once.
%
%==========================================================================
%==========================================================================

% NOTE: Numbered equations in the code (e.g., Eq. 21) correspond to 
%  equations in the [2] reference.

%==========================================================================
% Load GPS Accepted WGS-84 Constants 
%==========================================================================
muE = 3.986005e14;     % WGS-84 value, m^3/s^2
wE  = 7.2921151467e-5; % WGS-84 value, rad/s 
c   = 2.99792458e8;    % GPS acceptd speed of light, m/s
%PI = 3.1415926535898; % accepted GPS value for pi (not needed here)

%==========================================================================
% Initialize Output Variables for Speed 
%==========================================================================
sz         = size(t_input,1);
x          = ones(sz,3) * NaN;
v          = ones(sz,3) * NaN; 
health     = ones(sz,1) * NaN; 
satClkCorr = ones(sz,1) * NaN;
relcorr    = ones(sz,1) * NaN;

%==========================================================================
% Pull Out Correct Ephemerides 
%==========================================================================

% Pull out ephemerides for PRN in question
kk  = find(ephem_all(:,1) == prn);  % kk is vector containing row numbers of ephem_all that are for sat.no. 'index' 
sat_ephem = ephem_all(kk,:);        % sat_ephem is matrix of all ephem data for each entry of sat.no. 'index'

% No matching PRN found, returning data will be NaNs
if isempty(kk),return,end 

% Remove bad ephemerides. Occasionally, new ephemerides are uploaded at odd
%  times which occur prior to the 00:00 epochs and mostly commonly occur
%  at 59:44. Therefore, if these are present, we want to use the new info
%  and discard the ephemeris following it at 00:00.
% % % newuploads = find( mod( sat_ephem(:,17), 3600 ) ~= 0 );
% % % 
% % % if ~isempty(newuploads) % If new uploads exist
% % %     badrows = [];
% % %     for zz = 1:length(newuploads)  % Loop over new uploads
% % %         if (newuploads(zz) == length(sat_ephem(:,1))) % If new upload is last ephem entry
% % %             continue;
% % %         else
% % %             % Check that new upload is within 4 minutes of next ephem entry
% % %             %  In which case, remove the usual 00:00 entry and use the 
% % %             %  new upload
% % %             if ((sat_ephem(newuploads(zz)+1,17) - sat_ephem(newuploads(zz),17)) < 240)
% % %                 badrows = [badrows; newuploads(zz) + 1];    
% % %             end
% % %         end
% % %     end % END of FOR loop
% % %     % Remove unwanted ephemeris entries
% % %     sat_ephem(badrows,:) = [];
% % % end % END of IF new uploads exist

%==========================================================================
% Start Main Calculation Loop 
%==========================================================================
if strcmp(input,'eph')
    % Compute elapsed times of each ephemeris epoch wrt first entry, seconds
    dt_ephem = (sat_ephem(:,19) - sat_ephem(1,19))*604800 + (sat_ephem(:,17) - sat_ephem(1,17));


    % Compute elapsed times of each input time wrt first ephemeris entry, seconds
    dt_input = (t_input(:,1) - sat_ephem(1,19))*604800 + (t_input(:,2) - sat_ephem(1,17));
else
    % Modify t_sec so that they are all referenced to the gps week of ephemeris
    gpswk_ephem = ephem_all(1,19);  % GPS Week in YUMA almanac ephemeris matrix
end

for tt = 1:sz % loop through all input times

    if strcmp(input,'eph') 
    % Pull out most recent ephemeris values
        jj = max( find(dt_input(tt) >= dt_ephem) ); % sat_ephem(:,17) = toe (sec into GPS week) of each entry
        dt  = dt_input(tt) - dt_ephem(jj); % seconds difference from epoch                                      % jj = row of specific sat. ephem. data with epoch closest to input time
    elseif strcmp(input,'alm') 
        % for almanac there is only one value
        jj=1;
        t_sec = (t_input(tt,1) - gpswk_ephem).*7.*86400 + t_input(tt,2);
        toe = sat_ephem(jj,17);         % time of ephemeris
        dt  = t_sec - toe;               % seconds difference from epoch
        
    end
    % Pull out nearest ephemeris values                                                                                        
    %[mn,jj] = min(abs( dt_input(tt) - dt_ephem ));                                             
                                                      
    if isempty(jj),continue,end  % no matching ephemeris time found. continue to next input time 

    % Pull out common variables from the ephemeris matrix
    %======================================================================
    %toe = sat_ephem(jj,17);           % time of ephemeris
   
    
    a   = sat_ephem(jj,5)^2;           % semimajor axis, sqrt(a) = gps_ephem_all(:,5) (meters)
    ecc = sat_ephem(jj,4);             % eccentricity
    n0  = sqrt(muE/a^3);               % nominal mean motion (rad/s)
    n   = n0 + sat_ephem(jj,3);        % corrected mean motion, delta_n = gps_ephem_all(:,3)
    M   = sat_ephem(jj,2) + n*dt;      % mean anomaly, M0 = gps_ephem_all(:,2)

    % GPS Satellite Clock Correction, meters
    %======================================================================
    af0 = sat_ephem(jj,21);
    af1 = sat_ephem(jj,22);
    af2 = sat_ephem(jj,23);

    satClkCorr(tt,1) = c*((af2*dt + af1)*dt + af0); % meters

    % Compute perigee, true and eccentric anomaly...
    %======================================================================

    % Load argument of perigee to a local variable and add perigee rate, rad
    perigee  = sat_ephem(jj,8) + sat_ephem(jj,24) * dt;  

    % Compute Eccentric Anomaly, rad
    E    = mean2eccentric(M,ecc);
    cosE = cos(E);  
    sinE = sin(E);

    % Compute rate of change of Eccentric Anomaly, rad/s (Eq. 20)
    E_dot = n / (1-ecc*cosE);

    % Compute true anomaly, rad
    nu    = atan2( sqrt(1 - ecc*ecc).*sinE,  cosE-ecc ); 
    cosnu = cos(nu);  
    sinnu = sin(nu);  

    % Compute the argument of latitude, rad 
    phi = nu + perigee;  % true anomaly + argument of perigee

    % Compute corrections to argument of latitude, radius, and inclination
    %======================================================================
    costwophi = cos(2*phi);  
    sintwophi = sin(2*phi);

    delta_u = sat_ephem(jj,12) * sintwophi + ... % Cus = gps_ephem_all(jj,12)
              sat_ephem(jj,11) * costwophi;      % Cuc = gps_ephem_all(jj,11)

    delta_r = sat_ephem(jj,14) * sintwophi + ... % Crs = gps_ephem_all(jj,14)
              sat_ephem(jj,13) * costwophi;      % Crc = gps_ephem_all(jj,13)

    delta_i = sat_ephem(jj,16) * sintwophi + ... % Cis = gps_ephem_all(jj,16)
              sat_ephem(jj,15) * costwophi;      % Cic = gps_ephem_all(jj,15)

    u   = phi + delta_u;                                       % corrected argument of latitude
    r   = a * (1 - ecc*cosE) + delta_r;                        % corrected radius  
    inc = sat_ephem(jj,7) + delta_i + sat_ephem(jj,10) * dt;   % corrected inclination 
                                                               % i_dot = sat_ephem(jj,10)
    cosu = cos(u);  cos2u = cos(2*u);  
    sinu = sin(u);  sin2u = sin(2*u);

    % Compute Rates of Change of true anomaly, arg. of lat., radius, inclination
    %======================================================================

    % Compute rate of change of true anomaly, rad/s  
    %   used in   http://www.ngs.noaa.gov/gps-toolbox/bc_velo/bc_velo.c
    nu_dot = sinE*E_dot*(1+ecc*cosnu) / (sinnu*(1-ecc*cosE));

    % Eq 24 and used in skyplot.cpp 
    %nu_dot2 = a*a*sqrt(1 - e^2)*n ./ (a .* (1 - e .* cos(E))).^2; 

    % NOTE: the 2 previous equations for nu_dot are algebraically equal

    % Eq. 25, 26 and 24
    u_dot = nu_dot + 2*(sat_ephem(jj,12)*cos2u-sat_ephem(jj,11)*sin2u)*nu_dot;

    % Eq. 19, 20, 22 and 24 
    r_dot = a*ecc*sinE*n/(1-ecc*cosE) + 2*(sat_ephem(jj,14)*cos2u-sat_ephem(jj,13)*sin2u)*nu_dot;

    % Same format as Eq. 22 and 26 but with Cic and Cis instead
    i_dot = sat_ephem(jj,10) + 2*(sat_ephem(jj,16)*cos2u-sat_ephem(jj,15)*sin2u)*nu_dot;

    % Compute satellite position in orbital plane (Eq. 13)
    %======================================================================
    xo = r * cosu;    % satellite x-position in orbital plane
    yo = r * sinu;    % satellite y-position in orbital plane

    % Compute satellite velocity in orbital plane, Eq. 18
    %======================================================================
    xo_dot = r_dot*cosu - yo*u_dot;
    yo_dot = r_dot*sinu + xo*u_dot;

    % Corrected longitude of ascending node for node rate and Earth rotation
    %======================================================================
    % Ascending node = ephem_all(jj,6)
    node = sat_ephem(jj,6) + (sat_ephem(jj,9) - wE)*dt -  (wE * sat_ephem(jj,17)); % Toe = gps_ephem_all(jj,17)

    node_dot = sat_ephem(jj,9) - wE;    %Eq. 10,  node rate = ephem_all(jj,9)

    % Calculate GPS Satellite Position in ECEF (m)
    %======================================================================
    cosi = cos(inc);    sini = sin(inc);
    coso = cos(node);   sino = sin(node);


    % Satellite position in ECEF (m)
    x(tt,1) = xo*coso - yo*cosi*sino;  %x-position  

    x(tt,2) = xo*sino + yo*cosi*coso;  %y-position 

    x(tt,3) = yo*sini;                 %z-position


    % Calculate Satellite Velocity in ECEF (m/s)
    %======================================================================

    % Full velocity expression, Eq. 9
    %  Also presented in  http://www.ngs.noaa.gov/gps-toolbox/bc_velo/bc_velo.c
    v(tt,1) = (xo_dot - yo*cosi*node_dot)*coso - (xo*node_dot + yo_dot*cosi - yo*sini*i_dot)*sino;

    v(tt,2) = (xo_dot - yo*cosi*node_dot)*sino + (xo*node_dot + yo_dot*cosi - yo*sini*i_dot)*coso;

    v(tt,3) = yo_dot*sini + yo*cosi*i_dot;

    % Calculate relativistic correction (p. 93 of IS-GPS-200G)
    %======================================================================
    relcorr(tt,1) = c * (4.442807633e-10 * ecc * sat_ephem(jj,5) * sinE); % meters

    % Keep track of health of each satellite
    %======================================================================      
    health(tt,1) = sat_ephem(jj,25); % satellite health (0.00 is useable)

end % END of t_input loop =================================================
%==========================================================================    
end

function [E] = mean2eccentric(M,ecc) 

%==========================================================================
%==========================================================================
% [E] = mean2eccentric(M,e)
%
% Calculates the eccentric anomaly from the mean anomaly using Newton's 
%  method. The tolerance is set to 1e-12.
%
%
% Author: Ben K. Bradley
% Last Revision Date: 26 October 2010
%
%
%
% INPUT:            Description                                       Units
%
%  M          - mean anomaly (single value or vector)                   rad
%  ecc        - eccentricity orbit (single value or vector)      
%
% OUTPUT:
%
%  E          - eccentric anomaly (same size as input M)                rad
%
%
% Coupling:
%
%  none
%
% References:
%
%  [1] Curtis, H.D. "Orbital Mechanics for Engineering Students".
%          Elsevier Ltd. 2005.
%
%==========================================================================
%==========================================================================

if any(find((ecc > 0.999))) || any(find((ecc < 0)))
    errordlg({'Eccentric anomaly not solved for in mean2eccentric.m','Not an elliptic orbit.'},'Error!');
end

%==========================================================================

% Set tolerance
twopi   = 2*pi;
tol     = 1.0e-12;
maxiter = 20;


% Make sure mean anomaly is positive and within 2pi
M = rem(M,twopi);

M = (M<0)*twopi + M;


%  Make an initial guess for E, Eccentric Anomaly
% if (M < pi)
%     E = M + ecc/2;
% else
%     E = M - ecc/2;
% end
sinM = sin(M);

E = M + (ecc.*sinM) ./ (1 - sin(M + ecc) + sinM);


% Initialize iteration count and error
iter = 1;
err  = 1;


% Iterate to find E, eccentric anomaly
while any(find(abs(err) > tol)) && (iter <= maxiter)
    
    err  = (E - ecc.*sin(E) - M) ./ (1 - ecc.*cos(E));

    E    = E - err;
    
    iter = iter + 1;
    
    if (iter > maxiter)
        warning('Iterations maxed out in mean2eccentric.m'); %#ok<WNTAG>
    end
    
end
end