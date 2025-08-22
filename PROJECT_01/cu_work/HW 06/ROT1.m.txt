function [ c1 ] = ROT1( theta )
%This function preforms a rotation about the 1-axis given an input of theta
%in degrees

    c1 = [1, 0, 0;
          0, cosd(theta), sind(theta);
          0, -1*sind(theta), cosd(theta)];

end

