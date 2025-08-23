% ASEN 5090 HW # 5
% Author: CG

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
[brdc2630,ionoparams2630] = read_GPSbroadcast('brdc2630.16n');
[brdc2640,ionoparams2640] = read_GPSbroadcast('brdc2640.16n');
[brdc2650,ionoparams2650] = read_GPSbroadcast('brdc2650.16n');
[brdc2660,ionoparams2660] = read_GPSbroadcast('brdc2660.16n');
[brdc2670,ionoparams2670] = read_GPSbroadcast('brdc2670.16n');
[brdc2680,ionoparams2680] = read_GPSbroadcast('brdc2680.16n');

% yumaEpoch = [gps_week, gps_seconds]
yumaEpoch = [gps_yuma(1,19), gps_yuma(1,17)];
yumaEpochJD = gpsTime2JD(yumaEpoch(1), yumaEpoch(2));
yumaEpochJD = yumaEpochJD + 17/86400;

% 1. Take a look at the provided YUMA almanac file, what PRN's are in which
% orbit plane? Identify what orbit plane each PRN is in and add it to the
% table in your report.

    orbitPlane = rad2deg(gps_yuma(:,6));
    yumaGst0 = jd2gst(yumaEpochJD);
    for ii = 1:1:length(orbitPlane)
        orbitPlane(ii) = orbitPlane(ii) + yumaGst0;
        aa(ii) = orbitPlane(ii);
        if orbitPlane(ii) < 0
            orbitPlane(ii) = orbitPlane(ii) + 360;
        elseif orbitPlane(ii) > 360
            orbitPlane(ii) = orbitPlane(ii) - 360;
        end
    end

    disp('Solution to problem 1')
    orbitPlane
    
% 2. Write a function lla2ecef.m that converts from geodetic lat, lon and
% ellipsoid height to x, y, z in the Earth Centered Fixed coordinate frame.
% Make sure that it takes the WGS-84 ellipsoid parameters into account.
% Provide the ECEF coordinates for Boulder, CO (40, -105, 1631) in a table.
    LatBold = 40; LongBold = -105; AltBold = 1631;
    obsECEF = LLA2ECEF( LatBold, LongBold, AltBold, wgs84 )
    %obsECEF = lla2ecef([40, -105, 1631])
    
    disp('Solution to problem 2')
    obsECEF
        
% 3. Write a function ecef2azelrange.m that returns azimuth, elevation and
% range. Provide a table of all PRN's visible (above elevation 0) at 2 PM
% MST on September 18th, in Boulder, CO (40, -105, 1631), and their
% respective azimuths elevations and ranges with 2 digits of precision.
    lat = 40; long = -105;
    q3 = brdc2620(347:378,:);
    satJD = gpsTime2JD(1915, 72000)
    GST = jd2gst(satJD)
    for ii = 1:1:length(q3)
        M = q3(ii,2);
        ecc = q3(ii,4);
        a = (q3(ii,5))^2;
        a = a*1E-3;
        LOA = rad2deg(q3(ii,6));
        Om = LOA + GST;      
        om = rad2deg(q3(ii,8));
        E = mean2eccentric(M,ecc);
        tanNuo2 = sqrt((1 + ecc)/(1 - ecc))*tan(E/2);       
        Nuo2 = atan(tanNuo2);                           
        
        nu = rad2deg(2*Nuo2);                            % degrees
        if (nu < 0)
            nu = nu + 360;                               % quadrant check
        end
        
        i = rad2deg(q3(ii,7));
        
        [ R_ijk, V_ijk ] = COE2RV( a, ecc, i, Om, om, nu, wgs84.GM );
        
        satECEF = ROT3(GST)*R_ijk;
        rad(ii,:) = R_ijk;
        norm(R_ijk,2);
        
        [range, az, ele] = ecef2topo( obsECEF, satECEF, lat, long, wgs84 );
        prnAzElRa(ii,:) = [ii, az, ele, range];
    end
    
    disp('Solution to problem 3')
    prnAzElRa
    
% 4. Generate azimuth-elevation plots for visible satellites (above 0
% degree elevation) for the entire 24 hours (GPS time) on September 18th
% 2016. Do this for the following locations:
%   -   Boulder, CO (40,-105,1631)
%   -    0 N,  0 E
%   -   90 N,  0 E
%   -   Describe and discuss the differences in visibility.

    lat = 40; long = -105;
    q4 = brdc2620;
    for ii = 1:1:length(q4)
        satJD = gpsTime2JD(q4(ii,19), q4(ii,20));
        GST = jd2gst(satJD);
        M = q4(ii,2);
        ecc = q4(ii,4);
        a = (q4(ii,5))^2;
        LOA = rad2deg(q4(ii,6));
        Om = LOA + GST;      
        om = rad2deg(q4(ii,8));
        E = mean2eccentric(M,ecc);
        tanNuo2 = sqrt((1 + ecc)/(1 - ecc))*tan(E/2);       
        Nuo2 = atan(tanNuo2);                           
        
        nu = rad2deg(2*Nuo2);                            % degrees
        if (nu < 0)
            nu = nu + 360;                               % quadrant check
        end
        
        i = rad2deg(q4(ii,7));
        
        [ R_ijk, V_ijk ] = COE2RV( a, ecc, i, Om, om, nu, wgs84.GM );
        
        satECEF = ROT3(GST)*R_ijk;
        rad(ii,:) = R_ijk;
        norm(R_ijk,2);
        
        [range, az, ele] = ecef2topo( obsECEF, satECEF, lat, long, wgs84 );
        prnAzElRa2(ii,:) = [q4(ii,1), az, ele, range];
    end
    check = prnAzElRa2;
    vv = length(prnAzElRa2);
for jj = 1:1:length(check)
    mm = check(jj,3);
    if mm > 0
        prnPlotter(jj,:) = prnAzElRa2(jj,:);
    end
end

sorPrnPlo = sortrows(prnPlotter,1);

p01 = sorPrnPlo(269:274,:);
p02 = sorPrnPlo(275:281,:);
p03 = sorPrnPlo(282:287,:);
p04 = sorPrnPlo(288:291,:);
p05 = sorPrnPlo(292:296,:);
p06 = sorPrnPlo(297:301,:);
p07 = sorPrnPlo(302:305,:);
p08 = sorPrnPlo(306:311,:);
p09 = sorPrnPlo(312:318,:);
p10 = sorPrnPlo(319:322,:);
p11 = sorPrnPlo(323:327,:);
p13 = sorPrnPlo(328:334,:);
p14 = sorPrnPlo(335:341,:);
p15 = sorPrnPlo(342:348,:);
p17 = sorPrnPlo(349:352,:);
p18 = sorPrnPlo(353:358,:);
p19 = sorPrnPlo(359:365,:);
p20 = sorPrnPlo(366:371,:);
p21 = sorPrnPlo(372:375,:);
p22 = sorPrnPlo(376:381,:);
p23 = sorPrnPlo(382:386,:);
p24 = sorPrnPlo(387:393,:);
p27 = sorPrnPlo(394:397,:);
p29 = sorPrnPlo(398:401,:);
p30 = sorPrnPlo(402:405,:);
p31 = sorPrnPlo(406:409,:);
p32 = sorPrnPlo(410:416,:);


figure()
polarplot(deg2rad(p01(:,2)),(90 - p01(:,3)),'*b')
title('Sats over Boulder')
hold on
polarplot(deg2rad(p02(:,2)),(90 - p02(:,3)),'*b')
polarplot(deg2rad(p03(:,2)),(90 - p03(:,3)),'*b')
polarplot(deg2rad(p04(:,2)),(90 - p04(:,3)),'*b')
polarplot(deg2rad(p05(:,2)),(90 - p05(:,3)),'*b')
polarplot(deg2rad(p06(:,2)),(90 - p06(:,3)),'*b')
polarplot(deg2rad(p07(:,2)),(90 - p07(:,3)),'*b')
polarplot(deg2rad(p08(:,2)),(90 - p08(:,3)),'*b')
polarplot(deg2rad(p09(:,2)),(90 - p09(:,3)),'*b')
polarplot(deg2rad(p10(:,2)),(90 - p10(:,3)),'*b')
polarplot(deg2rad(p11(:,2)),(90 - p11(:,3)),'*b')
polarplot(deg2rad(p13(:,2)),(90 - p13(:,3)),'*b')
polarplot(deg2rad(p14(:,2)),(90 - p14(:,3)),'*b')
polarplot(deg2rad(p15(:,2)),(90 - p15(:,3)),'*b')
polarplot(deg2rad(p17(:,2)),(90 - p17(:,3)),'*b')
polarplot(deg2rad(p18(:,2)),(90 - p18(:,3)),'*b')
polarplot(deg2rad(p19(:,2)),(90 - p19(:,3)),'*b')
polarplot(deg2rad(p20(:,2)),(90 - p20(:,3)),'*b')
polarplot(deg2rad(p21(:,2)),(90 - p21(:,3)),'*b')
polarplot(deg2rad(p22(:,2)),(90 - p22(:,3)),'*b')
polarplot(deg2rad(p23(:,2)),(90 - p23(:,3)),'*b')
polarplot(deg2rad(p24(:,2)),(90 - p24(:,3)),'*b')
polarplot(deg2rad(p27(:,2)),(90 - p27(:,3)),'*b')
polarplot(deg2rad(p29(:,2)),(90 - p29(:,3)),'*b')
polarplot(deg2rad(p30(:,2)),(90 - p30(:,3)),'*b')
polarplot(deg2rad(p31(:,2)),(90 - p31(:,3)),'*b')
polarplot(deg2rad(p32(:,2)),(90 - p32(:,3)),'*b')

lat = 0.0; long = 0.0;
    q4 = brdc2620;
    for ii = 1:1:length(q4)
        satJD = gpsTime2JD(q4(ii,19), q4(ii,20));
        GST = jd2gst(satJD);
        M = q4(ii,2);
        ecc = q4(ii,4);
        a = (q4(ii,5))^2;
        LOA = rad2deg(q4(ii,6));
        Om = LOA + GST;      
        om = rad2deg(q4(ii,8));
        E = mean2eccentric(M,ecc);
        tanNuo2 = sqrt((1 + ecc)/(1 - ecc))*tan(E/2);       
        Nuo2 = atan(tanNuo2);                           
        
        nu = rad2deg(2*Nuo2);                            % degrees
        if (nu < 0)
            nu = nu + 360;                               % quadrant check
        end
        
        i = rad2deg(q4(ii,7));
        
        [ R_ijk, V_ijk ] = COE2RV( a, ecc, i, Om, om, nu, wgs84.GM );
        
        satECEF = ROT3(GST)*R_ijk;
        rad(ii,:) = R_ijk;
        norm(R_ijk,2);
        
        [range, az, ele] = ecef2topo( obsECEF, satECEF, lat, long, wgs84 );
        prnAzElRa3(ii,:) = [q4(ii,1), az, ele, range];
    end
    check2 = prnAzElRa3;
    vv2 = length(prnAzElRa3);
for jj = 1:1:length(check2)
    mm = check2(jj,3);
    if mm > 0
        prnPlotter2(jj,:) = prnAzElRa3(jj,:);
    end
end

sorPrnPlo2 = sortrows(prnPlotter2,1);

p01 = sorPrnPlo2(200:208,:);
p02 = sorPrnPlo2(209:214,:);
p03 = sorPrnPlo2(215:221,:);
p04 = sorPrnPlo2(222:227,:);
p05 = sorPrnPlo2(228:233,:);
p06 = sorPrnPlo2(234:240,:);
p07 = sorPrnPlo2(241:247,:);
p08 = sorPrnPlo2(248:254,:);
p09 = sorPrnPlo2(255:261,:);
p10 = sorPrnPlo2(262:267,:);
p11 = sorPrnPlo2(268:273,:);
p12 = sorPrnPlo2(274:280,:);
p13 = sorPrnPlo2(281:286,:);
p14 = sorPrnPlo2(287:292,:);
p15 = sorPrnPlo2(293:298,:);
p16 = sorPrnPlo2(299:305,:);
p17 = sorPrnPlo2(306:312,:);
p18 = sorPrnPlo2(313:319,:);
p19 = sorPrnPlo2(320:324,:);
p20 = sorPrnPlo2(324:331,:);
p21 = sorPrnPlo2(332:338,:);
p22 = sorPrnPlo2(339:345,:);
p23 = sorPrnPlo2(346:351,:);
p24 = sorPrnPlo2(352:361,:);
p25 = sorPrnPlo2(362:368,:);
p26 = sorPrnPlo2(369:375,:);
p27 = sorPrnPlo2(376:382,:);
p28 = sorPrnPlo2(383:388,:);
p29 = sorPrnPlo2(389:397,:);
p30 = sorPrnPlo2(398:403,:);
p31 = sorPrnPlo2(404:409,:);
p32 = sorPrnPlo2(410:416,:);

figure()
polarplot(deg2rad(p01(:,2)),(90 - p01(:,3)),'*b')
title('Sats over 0N 0E')
hold on
polarplot(deg2rad(p02(:,2)),(90 - p02(:,3)),'*b')
polarplot(deg2rad(p03(:,2)),(90 - p03(:,3)),'*b')
polarplot(deg2rad(p04(:,2)),(90 - p04(:,3)),'*b')
polarplot(deg2rad(p05(:,2)),(90 - p05(:,3)),'*b')
polarplot(deg2rad(p06(:,2)),(90 - p06(:,3)),'*b')
polarplot(deg2rad(p07(:,2)),(90 - p07(:,3)),'*b')
polarplot(deg2rad(p08(:,2)),(90 - p08(:,3)),'*b')
polarplot(deg2rad(p09(:,2)),(90 - p09(:,3)),'*b')
polarplot(deg2rad(p10(:,2)),(90 - p10(:,3)),'*b')
polarplot(deg2rad(p11(:,2)),(90 - p11(:,3)),'*b')
polarplot(deg2rad(p12(:,2)),(90 - p12(:,3)),'*b')
polarplot(deg2rad(p13(:,2)),(90 - p13(:,3)),'*b')
polarplot(deg2rad(p14(:,2)),(90 - p14(:,3)),'*b')
polarplot(deg2rad(p15(:,2)),(90 - p15(:,3)),'*b')
polarplot(deg2rad(p16(:,2)),(90 - p16(:,3)),'*b')
polarplot(deg2rad(p17(:,2)),(90 - p17(:,3)),'*b')
polarplot(deg2rad(p18(:,2)),(90 - p18(:,3)),'*b')
polarplot(deg2rad(p19(:,2)),(90 - p19(:,3)),'*b')
polarplot(deg2rad(p20(:,2)),(90 - p20(:,3)),'*b')
polarplot(deg2rad(p21(:,2)),(90 - p21(:,3)),'*b')
polarplot(deg2rad(p22(:,2)),(90 - p22(:,3)),'*b')
polarplot(deg2rad(p23(:,2)),(90 - p23(:,3)),'*b')
polarplot(deg2rad(p24(:,2)),(90 - p24(:,3)),'*b')
polarplot(deg2rad(p25(:,2)),(90 - p25(:,3)),'*b')
polarplot(deg2rad(p26(:,2)),(90 - p26(:,3)),'*b')
polarplot(deg2rad(p27(:,2)),(90 - p27(:,3)),'*b')
polarplot(deg2rad(p28(:,2)),(90 - p28(:,3)),'*b')
polarplot(deg2rad(p29(:,2)),(90 - p29(:,3)),'*b')
polarplot(deg2rad(p30(:,2)),(90 - p30(:,3)),'*b')
polarplot(deg2rad(p31(:,2)),(90 - p31(:,3)),'*b')
polarplot(deg2rad(p32(:,2)),(90 - p32(:,3)),'*b')


lat = 90.0; long = 0.0;
    q4 = brdc2620;
    for ii = 1:1:length(q4)
        satJD = gpsTime2JD(q4(ii,19), q4(ii,20));
        GST = jd2gst(satJD);
        M = q4(ii,2);
        ecc = q4(ii,4);
        a = (q4(ii,5))^2;
        LOA = rad2deg(q4(ii,6));
        Om = LOA + GST;      
        om = rad2deg(q4(ii,8));
        E = mean2eccentric(M,ecc);
        tanNuo2 = sqrt((1 + ecc)/(1 - ecc))*tan(E/2);       
        Nuo2 = atan(tanNuo2);                           
        
        nu = rad2deg(2*Nuo2);                            % degrees
        if (nu < 0)
            nu = nu + 360;                               % quadrant check
        end
        
        i = rad2deg(q4(ii,7));
        
        [ R_ijk, V_ijk ] = COE2RV( a, ecc, i, Om, om, nu, wgs84.GM );
        
        satECEF = ROT3(GST)*R_ijk;
        rad(ii,:) = R_ijk;
        norm(R_ijk,2);
        
        [range, az, ele] = ecef2topo( obsECEF, satECEF, lat, long, wgs84 );
        prnAzElRa4(ii,:) = [q4(ii,1), az, ele, range];
    end
    check3 = prnAzElRa4;
    vv3 = length(prnAzElRa4);
for jj = 1:1:length(check3)
    mm = check3(jj,3);
    if mm > 0
        prnPlotter4(jj,:) = prnAzElRa4(jj,:);
    end
end

sorPrnPlo4 = sortrows(prnPlotter4,1);

p01 = sorPrnPlo4(230:235,:);
p02 = sorPrnPlo4(236:242,:);
p03 = sorPrnPlo4(243:248,:);
p04 = sorPrnPlo4(249:254,:);
p05 = sorPrnPlo4(255:259,:);
p06 = sorPrnPlo4(260:266,:);
p07 = sorPrnPlo4(267:273,:);
p08 = sorPrnPlo4(274:280,:);
p09 = sorPrnPlo4(281:287,:);
p10 = sorPrnPlo4(288:291,:);
p11 = sorPrnPlo4(292:298,:);
p12 = sorPrnPlo4(299:304,:);
p13 = sorPrnPlo4(305:308,:);
p14 = sorPrnPlo4(309:315,:);
p15 = sorPrnPlo4(316:322,:);
p16 = sorPrnPlo4(323:327,:);
p17 = sorPrnPlo4(328:331,:);
p18 = sorPrnPlo4(332:337,:);
p19 = sorPrnPlo4(338:343,:);
p20 = sorPrnPlo4(344:349,:);
p21 = sorPrnPlo4(350:355,:);
p22 = sorPrnPlo4(356:361,:);
p23 = sorPrnPlo4(362:365,:);
p24 = sorPrnPlo4(366:373,:);
p25 = sorPrnPlo4(374:379,:);
p26 = sorPrnPlo4(380:386,:);
p27 = sorPrnPlo4(387:390,:);
p28 = sorPrnPlo4(391:395,:);
p29 = sorPrnPlo4(396:400,:);
p30 = sorPrnPlo4(401:405,:);
p31 = sorPrnPlo4(406:411,:);
p32 = sorPrnPlo4(412:416,:);

figure()
polarplot(deg2rad(p01(:,2)),(90 - p01(:,3)),'*b')
title('Sats over 90N 0E')
hold on
polarplot(deg2rad(p02(:,2)),(90 - p02(:,3)),'*b')
polarplot(deg2rad(p03(:,2)),(90 - p03(:,3)),'*b')
polarplot(deg2rad(p04(:,2)),(90 - p04(:,3)),'*b')
polarplot(deg2rad(p05(:,2)),(90 - p05(:,3)),'*b')
polarplot(deg2rad(p06(:,2)),(90 - p06(:,3)),'*b')
polarplot(deg2rad(p07(:,2)),(90 - p07(:,3)),'*b')
polarplot(deg2rad(p08(:,2)),(90 - p08(:,3)),'*b')
polarplot(deg2rad(p09(:,2)),(90 - p09(:,3)),'*b')
polarplot(deg2rad(p10(:,2)),(90 - p10(:,3)),'*b')
polarplot(deg2rad(p11(:,2)),(90 - p11(:,3)),'*b')
polarplot(deg2rad(p12(:,2)),(90 - p12(:,3)),'*b')
polarplot(deg2rad(p13(:,2)),(90 - p13(:,3)),'*b')
polarplot(deg2rad(p14(:,2)),(90 - p14(:,3)),'*b')
polarplot(deg2rad(p15(:,2)),(90 - p15(:,3)),'*b')
polarplot(deg2rad(p16(:,2)),(90 - p16(:,3)),'*b')
polarplot(deg2rad(p17(:,2)),(90 - p17(:,3)),'*b')
polarplot(deg2rad(p18(:,2)),(90 - p18(:,3)),'*b')
polarplot(deg2rad(p19(:,2)),(90 - p19(:,3)),'*b')
polarplot(deg2rad(p20(:,2)),(90 - p20(:,3)),'*b')
polarplot(deg2rad(p21(:,2)),(90 - p21(:,3)),'*b')
polarplot(deg2rad(p22(:,2)),(90 - p22(:,3)),'*b')
polarplot(deg2rad(p23(:,2)),(90 - p23(:,3)),'*b')
polarplot(deg2rad(p24(:,2)),(90 - p24(:,3)),'*b')
polarplot(deg2rad(p25(:,2)),(90 - p25(:,3)),'*b')
polarplot(deg2rad(p26(:,2)),(90 - p26(:,3)),'*b')
polarplot(deg2rad(p27(:,2)),(90 - p27(:,3)),'*b')
polarplot(deg2rad(p28(:,2)),(90 - p28(:,3)),'*b')
polarplot(deg2rad(p29(:,2)),(90 - p29(:,3)),'*b')
polarplot(deg2rad(p30(:,2)),(90 - p30(:,3)),'*b')
polarplot(deg2rad(p31(:,2)),(90 - p31(:,3)),'*b')
polarplot(deg2rad(p32(:,2)),(90 - p32(:,3)),'*b')


% 5. Show analytically what the highest elevation should be for a user at
% the Noth Pole.

    %%%%% This was completed by hand %%%%%

% 6. Provide a plot of the number of visible satellites (elevation > 10
% degrees) over Boulder, CO for 30 hours on the 18th of September 2016 (GPS
% time).

   lat = 40; long = -105;
    q6 = brdc2620;
    for ii = 1:1:length(q6)
        satJD = gpsTime2JD(q6(ii,19), q6(ii,20));
        GST = jd2gst(satJD);
        M = q6(ii,2);
        ecc = q6(ii,4);
        a = (q6(ii,5))^2;
        LOA = rad2deg(q6(ii,6));
        Om = LOA + GST;      
        om = rad2deg(q6(ii,8));
        E = mean2eccentric(M,ecc);
        tanNuo2 = sqrt((1 + ecc)/(1 - ecc))*tan(E/2);       
        Nuo2 = atan(tanNuo2);                           
        
        nu = rad2deg(2*Nuo2);                            % degrees
        if (nu < 0)
            nu = nu + 360;                               % quadrant check
        end
        
        i = rad2deg(q6(ii,7));
        
        [ R_ijk, V_ijk ] = COE2RV( a, ecc, i, Om, om, nu, wgs84.GM );
        
        satECEF = ROT3(GST)*R_ijk;
        rad(ii,:) = R_ijk;
        norm(R_ijk,2);
        
        [range, az, ele] = ecef2topo( obsECEF, satECEF, lat, long, wgs84 );
        prnAzElRa5(ii,:) = [q6(ii,1), az, ele, range];
    end
    
    check5 = prnAzElRa5;
    vv = length(prnAzElRa5);
for jj = 1:1:length(check5)
    mm = check5(jj,3);
    if mm > 10
        prnPlotter5(jj,:) = prnAzElRa5(jj,:);
    end
end

lat = 40; long = -105;
    q6b = brdc2630;
    for ii = 1:1:length(q6b)
        satJD = gpsTime2JD(q6b(ii,19), q6b(ii,20));
        GST = jd2gst(satJD);
        M = q6b(ii,2);
        ecc = q6b(ii,4);
        a = (q6b(ii,5))^2;
        LOA = rad2deg(q6b(ii,6));
        Om = LOA + GST;      
        om = rad2deg(q6b(ii,8));
        E = mean2eccentric(M,ecc);
        tanNuo2 = sqrt((1 + ecc)/(1 - ecc))*tan(E/2);       
        Nuo2 = atan(tanNuo2);                           
        
        nu = rad2deg(2*Nuo2);                            % degrees
        if (nu < 0)
            nu = nu + 360;                               % quadrant check
        end
        
        i = rad2deg(q6b(ii,7));
        
        [ R_ijk, V_ijk ] = COE2RV( a, ecc, i, Om, om, nu, wgs84.GM );
        
        satECEF = ROT3(GST)*R_ijk;
        rad(ii,:) = R_ijk;
        norm(R_ijk,2);
        
        [range, az, ele] = ecef2topo( obsECEF, satECEF, lat, long, wgs84 );
        prnAzElRa6(ii,:) = [q6b(ii,1), az, ele, range];
    end
    
    check6 = prnAzElRa6;
    vv = length(prnAzElRa6);
for jj = 1:1:length(check6)
    mm = check6(jj,3);
    if mm > 10
        prnPlotter6(jj,:) = prnAzElRa6(jj,:);
    end
end

y = [8,8,11,11,11,11,8,8,11,11,9,9,9,9,10,10,11,11,8,8,12,12,10,10,10,10,11,11,11,11];
x = [0:1:29];
figure()
plot(x,y)
title('Visible satellites over Boulder')
axis([0 29 0 12.5])
xlabel('Time [hours]')
ylabel('Number of Satellites visible')

% 7. The second step of this homework is to compute the geometric range.
% Write a function compute_range.m that takes as input the GPS ephemeris
% data, a PRN, the recieve time in GPS seconds of the week, and the assumed
% receiver position coordinates in WGS-84, and returns the expected
% geometric range with (range1) and without (range0) the time of flight
% correction (in meters).
%   [range0, range1] = compute_range(eph, PRN, t, userpos);
% Provide range 1 and range 0 for PRN 2, the 18th of October, 2 PM MST in a
% table.
    
    eph = brdc2620(348,:); prn = 2; t = [1915 72000];
    userpos = LLA2ECEF( 40.0, -105.0, 1631, wgs84 );
    
    [range0, range1] = compute_range(eph, prn, t, userpos, wgs84);
    
    disp('Solution to problem 7')
    range0
    range1

% 8. Compute range1, and range0 for PRN 2, using GPS ephemeris for the entire
% week (7 days starting at 18th of October, GPS time). For the following
% plots use range1 from ephemeris as your truth. Plot range1 together with
% range0 (not corrected for transmission time), as well as the original
% range from your ecef2azelrange function. How large is the correction from
% satellite clock error, relativity and time of flight?

% Use your compute_range function to calcualte range1 with GPS YUMA almanac
% as input, do this for the entire week and compare to the truth.

% Compute range1 using GPS ephemeris, but do not incorporate the 2 hour
% updates (just use the first ephemeris and propagate it forward in time),
% do this for the entire week and compare to the truth.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% notes for class Monday
% broadcast2vx[gpsWeek, second of week]

% 


