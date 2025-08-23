function [ range0, range1 ] = compute_range(eph, prn, t, userpos, p)
    
        tol = 1.0E-8;

        satJD = gpsTime2JD(t(1), t(2));
        GST = jd2gst(satJD);
        M = eph(2);
        ecc = eph(4);
        a = (eph(5))^2;
        LOA = rad2deg(eph(6));
        Om = LOA + GST;      
        om = rad2deg(eph(8));
        E = mean2eccentric(M,ecc);
        tanNuo2 = sqrt((1 + ecc)/(1 - ecc))*tan(E/2);       
        Nuo2 = atan(tanNuo2);                           
        
        nu = rad2deg(2*Nuo2);                            % degrees
        if (nu < 0)
            nu = nu + 360;                               % quadrant check
        end
        
        i = rad2deg(eph(7));
        
        [ R_ijk, V_ijk ] = COE2RV( a, ecc, i, Om, om, nu, p.GM );
        
        satECEF_Tr = ROT3(GST)*R_ijk;
        
        R = abs(satECEF_Tr - userpos);
        
        check = norm(R,2);
     
        range1 = R;
        
        Tr = t(2);
        
        k = 0.0;
        
  while check > tol      
        
        Tt = Tr - check/p.spdLite;
        
        satJD = gpsTime2JD(t(1), Tt);
        GST = jd2gst(satJD);
        M = eph(2);
        ecc = eph(4);
        a = (eph(5))^2;
        LOA = rad2deg(eph(6));
        Om = LOA + GST;      
        om = rad2deg(eph(8));
        E = mean2eccentric(M,ecc);
        tanNuo2 = sqrt((1 + ecc)/(1 - ecc))*tan(E/2);       
        Nuo2 = atan(tanNuo2);                           
        
        nu = rad2deg(2*Nuo2);                            % degrees
        if (nu < 0)
            nu = nu + 360;                               % quadrant check
        end
        
        i = rad2deg(eph(7));
        
        [ R_ijk, V_ijk ] = COE2RV( a, ecc, i, Om, om, nu, p.GM );
        
        satECEF_Tt = ROT3(GST)*R_ijk;
        
        R = abs(satECEF_Tt - userpos);
        
        check = abs(norm(R,2) - norm(range1,2));
        
        k = k+1;
  end
  
  range0 = R;

end

