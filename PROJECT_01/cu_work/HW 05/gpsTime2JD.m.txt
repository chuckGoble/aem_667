function [ JulianDate ] = gpsTime2JD( gpsWeek, gpsSecond )
   
    gpsEpochJD = 2444244.5;
    
    a = gpsWeek * 7;
    
    a = a + gpsEpochJD;
    
    % b = gpsWeek + a;
    
    c = gpsSecond/86400;
    
    JulianDate = a + c;


end

