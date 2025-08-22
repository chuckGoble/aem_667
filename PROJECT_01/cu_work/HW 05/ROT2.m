function [ c2 ] = ROT2( theta )
%This function preforms a rotation about the 2-Axis given an input of theta
%in degrees
    
    c2 = [cosd(theta), 0 -1*sind(theta);
          0, 1, 0;
          sind(theta), 0, cosd(theta)];


end

