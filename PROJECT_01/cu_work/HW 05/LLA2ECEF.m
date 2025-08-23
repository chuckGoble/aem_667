function [ rECEF ] = LLA2ECEF( Lat, Long, Alt, p )
% Transfrom geodetic LLA to ECEF
    
    radEar = p.ellipSemi;
    f = p.f;
    eccEarSq = 2*f - f^2;
    

    C = radEar/(sqrt(1-((eccEarSq)*((sind(Lat))^2))));
    S = (radEar*(1 - eccEarSq))/(sqrt(1-((eccEarSq)*((sind(Lat))^2))));
    
    rECEF = [(C + Alt)*cosd(Lat)*cosd(Long);
             (C + Alt)*cosd(Lat)*sind(Long);
             (S + Alt)*sind(Lat)];
         
end

