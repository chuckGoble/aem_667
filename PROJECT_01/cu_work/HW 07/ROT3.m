function [ c3 ] = ROT3( theta )
%This function preforms a rotation about the 3-Axis given an input of theta
%in degrees

    c3 = [cosd(theta), sind(theta), 0;
            -1*sind(theta), cosd(theta),0;
            0, 0, 1];
end

