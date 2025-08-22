% Fix HW #5

close all; clear all; clc;

close all; clear all; clc; format longeng

wgs84.ellipSemi = 6378137.0;
wgs84.reciFlat = 298.257223563;
wgs84.f = 1/wgs84.reciFlat;
wgs84.earAngVel = 7292115.0E-11;
wgs84.GM = 3986004.418E8;
wgs84.spdLite = 2.99792458E8;

meanSolarDay = 86400;
meanSiderealDay = 86164.09954;
gps_yuma = read_GPSyuma('yuma0891.061440');
[brdc2620,ionoparams2620] = read_GPSbroadcast('brdc2620.16n');

% yumaEpoch = [gps_week, gps_seconds]
yumaEpoch = [gps_yuma(1,19), gps_yuma(1,17)];
yumaEpochJD = gpsTime2JD(yumaEpoch(1), yumaEpoch(2));
yumaEpochJD = yumaEpochJD + 17/86400;

% 3. Write a function ecef2azelrange.m that returns azimuth, elevation and
% range. Provide a table of all PRN's visible (above elevation 0) at 2 PM
% MST on September 18th, in Boulder, CO (40, -105, 1631), and their
% respective azimuths elevations and ranges with 2 digits of precision.
    lat = 40; long = -105;
    
    [health,x,v,relcorr,satClkCorr] = broadcast2xv(gps_yuma,[1915 72017],29, 'alm');
    obsECEF = lla2ecef([lat long 1631]);
    satECEF = x;
    [range, az, ele] = ecef2topo( obsECEF', satECEF', lat, long, wgs84 );
    
    
    
% 7. The second step of this homework is to compute the geometric range.
% Write a function compute_range.m that takes as input the GPS ephemeris
% data, a PRN, the recieve time in GPS seconds of the week, and the assumed
% receiver position coordinates in WGS-84, and returns the expected
% geometric range with (range1) and without (range0) the time of flight
% correction (in meters).
%   [range0, range1] = compute_range(eph, PRN, t, userpos);
% Provide range 1 and range 0 for PRN 2, the 18th of October, 2 PM MST in a
% table.

[ range0, range1] = compute_rangeRev02(brdc2620, 02, [1915 72017], obsECEF, wgs84)
    
    
    
    
    
    
    
    
    
    