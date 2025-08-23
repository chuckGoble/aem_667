function [ range, az, ele ] = ecef2topo( obsECEF, satECEF,lat,long, p )
%This function comptes the range, elevation and azimuth of a satellite
%given ECEF position and the latitude, longitude and altitude of the
%tracking station.
    
    rho = satECEF - obsECEF;
    
    QXx =[-sind(long), cosd(long), 0;
            -sind(lat)*cosd(long), -sind(lat)*sind(long), cosd(lat);
            cosd(lat)*cosd(long), cosd(lat)*sind(long), sind(lat)]; 
    
    rhoTran = QXx*rho;
    
    range = norm(rhoTran,2);
    
    unitRho = rhoTran./range;
    
    ele = asind(unitRho(3)/(norm(unitRho,2)));
    
    sinAz = unitRho(1)/cosd(ele);
    cosAz = unitRho(2)/cosd(ele);
    
    az = atan2d(sinAz,cosAz);
    if az < 0
        az = az + 360;
    end
    
end

