function [ gst0 ] = jd2gst( JD )
% Convert from Juilan Date to theta_GMST

    Tut1 = (JD - 2451545)/36525;

    thetaGMST0 = 67310.54841 + (3600*876600 + 8640184.812866)*Tut1 ...
        + 0.093104*(Tut1*Tut1) - 6.2E-6*(Tut1*Tut1*Tut1);
    
    if thetaGMST0 > 0
        while thetaGMST0 > 86400
            thetaGMST0 = thetaGMST0 - 86400;
        end
   else
        while thetaGMST0 < -86400
            thetaGMST0 = thetaGMST0 + 86400;
        end
   end
   
   gst0 = thetaGMST0/240;
   
   if gst0 < 0
       gst0 = gst0 + 360;
   end

end

