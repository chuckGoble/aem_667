function [ range0, range1,k,R0,R1 ] = compute_range(eph, prn, t, userpos, wgs84)
    
    tRec = t(2);

    [health,x,v,relcorr,satClkCorr] = broadcast2xv(eph,t,prn, 'eph');
    
    rho = x' - userpos;
    R0 = rho;
    range0 = norm(rho,2);
    
    c = wgs84.spdLite;
    
    old = range0;
    k = 0;
    err = 1;
    
    while(err > 1E-12) && k < 5
        
        tTrans = tRec - old/c;
        [health,x,v,relcorr,satClkCorr] = broadcast2xv(eph,[t(1) tTrans],prn, 'eph');
        phi = wgs84.earAngVel*(tRec - tTrans);
        phi = rad2deg(phi);
        Q = ROT3(phi);
        x = Q*x';
        geoR = x - userpos;
        range = norm(geoR,2);
        new = range;
        err = abs(old - new);
        old = new;
        k = k + 1;
    end
    R1 = geoR;
    range1 = new;
    
end