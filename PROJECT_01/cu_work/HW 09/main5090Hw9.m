% ASEN 5050 HW # 9
% Author: CG

close all; clear all; clc;

wgs84.ellipSemi = 6378137.0;
wgs84.reciFlat = 298.257223563;
wgs84.f = 1/wgs84.reciFlat;
wgs84.earAngVel = 7292115.0E-11;
wgs84.GM = 3986004.418E8;
wgs84.spdLite = 2.99792458E8;

meanSolarDay = 86400;
meanSiderealDay = 86164.09954;

% Set the trial number!
% Trial 1 = AMC station
% Trial 2 = Mystery station
trial = 1;

if trial == 1
    [fid, obsECEF, observables] = read_rinex_header('amc23200.16o');
    rinexv3 = read_rinex_obs5('amc23200.16o',1:1:32,86400);
    clc;
    data = rinexv3.data;
    val = (data(end,2)-data(1,2))/30;
    startTime = data(1,2) + 30;
else
    [fid, obsECEF, observables] = read_rinex_header('mystery3200.16o');
    rinexv3 = read_rinex_obs5('mystery3200.16o',1:1:32,86400);
    clc;
    data = rinexv3.data;
    startTime = data(1,2) + 60;
    val = (data(end,2)-startTime)/30;
end

dummy = 1;
storePreFitResVec = [];
storeDops = zeros(val,5);
storeDx = zeros(val,4);
storeDxENU = zeros(val,4);
storePostFitVec = [];
storePosENU = zeros(val,4);
storePosECEF = zeros(val,4);
storeDy = zeros(val,33);
storeEle = zeros(val,33);
storePostDy = zeros(val,33);

for jj = startTime:30:data(end,2)
    epoch = jj;
    aa = find(data(:,2) == epoch);
    newDat = data(aa,:);
    
    [brdc3000,ionoparams3000] = read_GPSbroadcast('brdc3200.16n');
        
        if (trial ~= 1) && (jj==startTime)
            p3 = ((1575.42^2)/(1575.42^2 - 1227.6^2))*newDat(:,7) -...
            ((1227.6^2)/(1575.42^2 - 1227.6^2))*newDat(:,8);
        end
    
    % Iterate for the first epoch of the mystery file
        if (trial ~=1)&&(jj==startTime)
            for kl = 1:1:10
                for ii = 1:1:size(newDat,1)
                    [health,x,v,relcorr,satClkCorr] = broadcast2xv(brdc3000,[data(1) epoch],newDat(ii,3), 'eph');data(1,2);
                    satECEF = x;
                    lla = ecef2lla(obsECEF');
                    [range, az,ele] = ecef2topo( obsECEF, satECEF', lla(1), lla(2), wgs84 );
                    [range0, range1,] = compute_range(brdc3000, newDat(ii,3), [data(1) epoch], obsECEF, wgs84);
                    rangeTopo(ii,:) = range;
                    azimuth(ii,:) = az;
                    elevation(ii,:) = ele;
                    geoRange(ii,:) = range0;
                    expRange(ii,:) = range1;
                    relativity(ii,:) = relcorr;
                    satClock(ii,:) = satClkCorr;
                    satPos(ii,:) = x;
                    clearvars satHealth
                    satHealth(ii,:) = health;
                end
                
                r1 = isnan(elevation);
                r2 = elevation < 10;
                r3 = satHealth ~= 0;
    
                tFall = r1 | r2 | r3;
                elevation(tFall,:) = [];
                geoRange(tFall,:) = [];
                expRange(tFall,:) = [];
                relativity(tFall,:) = [];
                satClock(tFall,:) = [];
                satPos(tFall,:) = [];
                rangeTopo(tFall,:) = [];
                azimuth(tFall,:) = [];
                newDat(tFall,:) = [];
                p3(tFall,:) = [];

                for kk = 1:1:size(elevation,1)
                    A(kk,:) = [-(satPos(kk,1)-obsECEF(1))/geoRange(kk), -(satPos(kk,2)-obsECEF(2))/geoRange(kk), -(satPos(kk,3)-obsECEF(3))/geoRange(kk), 1];
                end
            
                tropErr = 2.5./(sind(elevation));
                dy = p3 - expRange + satClock - relativity - tropErr;
            
                dx = inv(transpose(A)*A)*transpose(A)*dy;
                obsECEF = dx(1:3) + obsECEF;
                %ecef2lla(obsECEF')
                clear A dy
            end
        end
    
%     
        
    % Expected range in meters
    for ii = 1:1:size(newDat,1)
        [health,x,v,relcorr,satClkCorr] = broadcast2xv(brdc3000,[data(1) epoch],newDat(ii,3), 'eph');data(1,2);
        satECEF = x;
        lla = ecef2lla(obsECEF');
        [range, az,ele] = ecef2topo( obsECEF, satECEF', lla(1), lla(2), wgs84 );
        [range0, range1,] = compute_range(brdc3000, newDat(ii,3), [data(1) epoch], obsECEF, wgs84);
        rangeTopo(ii,:) = range;
        azimuth(ii,:) = az;
        elevation(ii,:) = ele;
        geoRange(ii,:) = range0;
        expRange(ii,:) = range1;
        relativity(ii,:) = relcorr;
        satClock(ii,:) = satClkCorr;
        satPos(ii,:) = x;
        satHealth(ii,:) = health;
        
    end     
        
    r1 = isnan(elevation);
    r2 = elevation < 10;
    r3 = satHealth ~= 0;
    
    tFall = r1 | r2 | r3;
    elevation(tFall,:) = [];
    geoRange(tFall,:) = [];
    expRange(tFall,:) = [];
    relativity(tFall,:) = [];
    satClock(tFall,:) = [];
    satPos(tFall,:) = [];
    rangeTopo(tFall,:) = [];
    azimuth(tFall,:) = [];
    newDat(tFall,:) = [];
    
    % Ionosphere free pseudorange observable P3
        p3 = ((1575.42^2)/(1575.42^2 - 1227.6^2))*newDat(:,6) -...
            ((1227.6^2)/(1575.42^2 - 1227.6^2))*newDat(:,7);
        
        if trial ~= 1
            p3 = ((1575.42^2)/(1575.42^2 - 1227.6^2))*newDat(:,7) -...
            ((1227.6^2)/(1575.42^2 - 1227.6^2))*newDat(:,8);
        end
    
    
    % Create the rotation matrix Q from ECEF to ENU
        qMat = ROT1(90 - lla(1))*ROT3(90+lla(2));
        qA = [qMat, [0;0;0]];
        qA = [qA; [0,0,0,1]];
    
    % Construct and print the A-Matrix as discussed in lecture or in your book.
        for kk = 1:1:size(elevation,1)
            A(kk,:) = [-(satPos(kk,1)-obsECEF(1))/geoRange(kk), -(satPos(kk,2)-obsECEF(2))/geoRange(kk), -(satPos(kk,3)-obsECEF(3))/geoRange(kk), 1];
        end
        % Rotate the A matrix from ECEF to ENU
            for ij = 1:1:size(A,1)
                temp = qA*A(ij,:)';
                A_ENU(ij,:) = temp;
            end
            
    % Compute and print the prefit residuals (dy) correcting for ionospheric errors, satellite clock, and relativity. 
        tropErr = 2.5./(sind(elevation));
        dy = p3 - expRange + satClock - relativity - tropErr;
        
    % Check dy for large number > 500
        r4 = abs(dy) > 500;
        dy(r4,:) = [];
        A(r4,:) = [];
        A_ENU(r4,:) = [];
        r5 = isnan(dy);
        dy(r5,:) = [];
        A(r5,:) = [];
        A_ENU(r5,:) = [];
        r6 = isnan(A);
        dy(r6,:) = [];
        A(r6,:) = [];
        A_ENU(r6,:) = [];
        r7 = isnan(A_ENU);
        dy(r7,:) = [];
        A(r7,:) = [];
        A_ENU(r7,:) = [];
    
    % Compute the least squares solution
        dx = inv(transpose(A)*A)*transpose(A)*dy;
        posECEF = dx(1:3) + obsECEF;
        dx_ENU = qA*dx;
        obsENU = qMat*obsECEF;
        posENU = dx_ENU(1:3) + obsENU;
    % Calculate the postfit residuals
        postFit = dy - A*dx;
        
    % Calculate the H matrix
        Hmat = inv(A_ENU'*A_ENU);
        dops = [Hmat(1,1),Hmat(2,2),Hmat(3,3),Hmat(4,4)];
    
    
        
        % Store all of the required data
        storePreFitResVec = [storePreFitResVec dy'];
        storeDx(dummy,:) = dx';
        storeDxENU(dummy,:)= dx_ENU';
        storePostFitVec = [storePostFitVec postFit'];
        storeDops(dummy,:) = [jj, dops, storeDops(dummy, length(dops)+2:end)];
        storePosECEF(dummy,:) = [jj,posECEF'];
        storePosENU(dummy,:) = [jj,posENU'];
        storeEle(dummy,:) = [elevation',storeEle(dummy, length(elevation)+1:end)];
        storeDy(dummy,:) =[dy', storeDy(dummy, length(dy) + 1:end)];
        storeDyPost(dummy,:) =[postFit', storeDy(dummy, length(postFit) + 1:end)];
        
        endTime = jj;
        
        clearvars -except storePreFitRes storeDx storePostFit dummy data obsECEF wgs84 storeDops trial startTime endTime storePosECEF storePosENU storeDy storeDyPost storePreFitResVec storePostFitVec dummy storeEle storeDxENU
        clearvars ii
        
        dummy = dummy + 1
   
end


    figure()
    subplot(3,1,1)
    % x-axis
    plot(storeDops(:,1)./3600,storeDxENU(:,1))
    title('East Relative Position verse Time')
    ylabel('Relative Position [m]')
    xlabel('Time [hr]')
    subplot(3,1,2)
    % y-axis
    plot(storeDops(:,1)./3600,storeDxENU(:,2))
    title('North Relative Position verse Time')
    ylabel('Relative Position [m]')
    xlabel('Time [hr]')
    subplot(3,1,3)
    % z-axis
    plot(storeDops(:,1)./3600,storeDxENU(:,3))
    title('Up Relative Position verse Time')
    ylabel('Relative Position [m]')
    xlabel('Time [hr]')
    
    figure()
    for ll=1:1:size(storeDy,1)
        plot(storeEle(ll,:),storeDy(ll,:),'b*')
        title('Prefit Residual versus Elevation Angle')
        xlabel('Elevation Angle [Deg]')
        ylabel('Prefit Residual [m]')
        hold on
    end
    
    figure()
    for ll=1:1:size(storeDy,1)
        plot(storeEle(ll,:),storeDyPost(ll,:),'b*')
        title('Postfit Residual versus Elevation Angle')
        xlabel('Elevation Angle [Deg]')
        ylabel('Prefit Residual [m]')
        hold on
    end
    figure()
    for ll=1:1:size(storeDy,1)
        plot(storeDops(ll,1)./3600,storeDyPost(ll,:),'b*')
        title('Postfit Residual versus Time [hr]')
        xlabel('Time [hr]')
        ylabel('Prefit Residual [m]')
        hold on
    end
    figure()
    for ll=1:1:size(storeDy,1)
        plot(storeDops(ll,1)./3600,storeDy(ll,:),'b*')
        title('Prefit Residual versus Time')
        xlabel('Time [hr]')
        ylabel('Prefit Residual [m]')
        hold on
    end
    
    figure()
    
    subplot(3,1,1)
    % East DOP
    plot(storeDops(:,1)./3600,storeDops(:,2))
    title('East DOP verse Time')
    ylabel('DOP [m]')
    xlabel('Time [hr]')
    hold on
    subplot(3,1,1)
    % Nort DOP
    plot(storeDops(:,1)./3600,storeDops(:,3))
    title('North DOP verse Time')
    ylabel('DOP [m]')
    xlabel('Time [hr]')
    subplot(3,1,1)
    % Up DOP
    plot(storeDops(:,1)./3600,storeDops(:,4))
    title('Up DOP verse Time')
    ylabel('DOP [m]')
    xlabel('Time [hr]')

% Find the rms for relative position
    rmsRelPos = rms(storeDx);
    
% Find the standard Deviation for the relative position
    stdRelX = std(storeDx(:,1));
    stdRelY = std(storeDx(:,2));
    stdRelZ = std(storeDx(:,3));
    
% Find the standard Deviation for the preFit Residuals
    stdPreFit = std(storePreFitResVec);
    stdPostFit = std(storePostFitVec);
    
% Find the RMS of the PreFit and PostFit Residuals
    rmsPrefit = rms(storePreFitResVec);
    rmsPostFit = rms(storePostFitVec);



