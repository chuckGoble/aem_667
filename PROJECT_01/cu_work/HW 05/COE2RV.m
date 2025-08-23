function [ R_ijk, V_ijk ] = COE2RV( a, e, i, Om, om, nu, mu )
% Author: Charles Goble

    p = a*(1-e^2);
    r = (p)/(1 + (e*cosd(nu)));
    
    Rp = [r*cosd(nu), r*sind(nu), 0]';
    Vp = (sqrt(mu/p))*[-1*sind(nu), (e + cosd(nu)), 0 ]';
    
    
    Q1 = [cosd(om), sind(om), 0; -1*sind(om), cosd(om), 0; 0, 0, 1];
    Q2 = [1, 0, 0; 0, cosd(i), sind(i); 0, -1*sind(i), cosd(i)];
    Q3 = [cosd(Om), sind(Om), 0; -1*sind(Om), cosd(Om), 0; 0, 0, 1];
    
    Q = Q1*Q2*Q3;
    Q = inv(Q);
    R_ijk = Q*Rp;
    
    V_ijk = Q*Vp;

end

