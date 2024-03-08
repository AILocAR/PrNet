function gpsPvt = GpsWlsPvtEKF(gnssMeas, allGpsEph, iono, Flag_discon, bRaw)
%gpsPvt = GpsWlsPvt(gnssMeas,allGpsEph,bRaw)
%compute PVT from gnssMeas
% Input: gnssMeas, structure of pseudoranges, etc. from ProcessGnssMeas
%        allGpsEph, structure with all ephemeris
%        [bRaw], default true, true => use raw pr, false => use smoothed
%
% Output:
% gpsPvt.FctSeconds    Nx1 time vector, same as gnssMeas.FctSeconds
%       .allLlaDegDegM Nx3 matrix, (i,:) = [lat (deg), lon (deg), alt (m)]
%       .sigmaLlaM     Nx3 standard deviation of [lat,lon,alt] (m)
%       .allBcMeters   Nx1 common bias computed with llaDegDegM
%       .allVelMps     Nx3 (i,:) = velocity in NED coords
%       .sigmaVelMps   Nx3 standard deviation of velocity (m/s)
%       .allBcDotMps   Nx1 common freq bias computed with velocity
%       .numSvs        Nx1 number of satellites used in corresponding llaDegDegM
%       .hdop          Nx1 hdop of corresponding fix
%
%Algorithm: Weighted Least Squares

%Author: Frank van Diggelen
%Modified by Xu Weng
%Open Source code for processing Android GNSS Measurements

if nargin < 5 % Number of function input arguments
    bRaw = true;
else
    %check that smoothed pr fields exists in input
    if any(~isfield(gnssMeas,{'PrSmM','PrSmSigmaM'}))
        error('If bRaw is false, gnssMeas must have fields gnssMeas.PrSmM and gnssMeas.PrSmSigmaM')
    end
end

xo =zeros(8,1);%initial state: [center of the Earth, bc=0, velocities = 0]'
xo(1:3) = [-1509706.868, 6195107.314,149731.63]';
% xo_MHE_k_N = xo;%initial state at the k-Nth epoch for MHE: [center of the Earth, bc=0, velocities = 0]'
xo_MHE = xo;%initial state for MHE: [center of the Earth, bc=0, velocities = 0]'

weekNum     = floor(gnssMeas.FctSeconds/GpsConstants.WEEKSEC);
%TBD check for week rollover here (it is checked in ProcessGnssMeas, but
%this function should stand alone, so we should check again, and adjust
%tRxSeconds by +- a week if necessary)
%btw, Q. why not just carry around fct and not worry about the hassle of
%weeknumber, and the associated week rollover problems?
% A. because you cannot get better than 1000ns (1 microsecond) precsision
% when you put fct into a double. And that would cause errors of ~800m/s * 1us
% (satellite range rate * time error) ~ 1mm in the range residual computation
% So what? well, if you start processing with carrier phase, these errors
% could accumulate.

N = length(gnssMeas.FctSeconds);
M = size(gnssMeas.PrM, 2);
gpsPvt.FctSeconds      = gnssMeas.FctSeconds;
gpsPvt.allXyzMMM       = zeros(N,3)+NaN;
gpsPvt.allLlaDegDegM   = zeros(N,3)+NaN;
gpsPvt.sigmaLLaM       = zeros(N,3)+NaN;
gpsPvt.allBcMeters     = zeros(N,1)+NaN;
gpsPvt.allVelMps       = zeros(N,3)+NaN;
gpsPvt.allVelMpsXyz    = zeros(N,3)+NaN;
gpsPvt.sigmaVelMps     = zeros(N,3)+NaN;
gpsPvt.allBcDotMps     = zeros(N,1)+NaN;
gpsPvt.numSvs          = zeros(N,1);
gpsPvt.hdop            = zeros(N,1)+inf;
gpsPvt.SvPosR          = [];
gpsPvt.SvPosRS         = [];
gpsPvt.Wpr             = cell(N,1);
gpsPvt.dtsv            = zeros(size(gnssMeas.PrM))+NaN;
gpsPvt.RrMps           = zeros(size(gnssMeas.PrM))+NaN;
gpsPvt.SvXyzMMM        = zeros(M,3,N)+NaN;
gpsPvt.iValid           = cell(N,1);

% Outputs of EKF
gpsPvt.allXp           = zeros(N,8)+NaN;
gpsPvt.allPp           = zeros(8,8,N)+NaN;
gpsPvt.allXhat         = zeros(N,8)+NaN;
gpsPvt.allPhat         = zeros(8,8,N)+NaN;
gpsPvt.allXhatS        = zeros(N,8)+NaN;
gpsPvt.allPhatS        = zeros(8,8,N)+NaN;
gpsPvt.allXyzMMMEKF    = zeros(N,3)+NaN;
gpsPvt.allLlaDegDegMEKF= zeros(N,3)+NaN;
gpsPvt.allBcMetersEKF  = zeros(N,1)+NaN;
gpsPvt.allfcMpsEKF     = zeros(N,1)+NaN;
gpsPvt.allVelMpsEKF    = zeros(N,3)+NaN;
gpsPvt.allXyzMMMEKFS    = zeros(N,3)+NaN;
gpsPvt.allLlaDegDegMEKFS= zeros(N,3)+NaN;
gpsPvt.allBcMetersEKFS  = zeros(N,1)+NaN;
gpsPvt.allfcMpsEKFS     = zeros(N,1)+NaN;
gpsPvt.allVelMpsEKFS    = zeros(N,3)+NaN;

% Outputs of MHE
gpsPvt.allXyzMMMMHE    = zeros(N,3)+NaN;
gpsPvt.allLlaDegDegMMHE= zeros(N,3)+NaN;
gpsPvt.allBcMetersMHE  = zeros(N,1)+NaN;
gpsPvt.allBcDotMpsMHE  = zeros(N,1)+NaN;
gpsPvt.allVelMpsMHE    = zeros(N,3)+NaN;

% Elapsed Time
elapsed_time_ekf = 0;
elapsed_time_rts = 0;

% GT_data = load('DenoisedPrM.csv');
% GT_data = load('PrM_Bias.csv');
% GT_data = load('PrM_Bias_2020-05-14-US-MTV-1_AR3_R.csv');
% GT_data = load('SvPVT3D_Error_label_dynamic_data_2020-06-05-US-MTV-1.csv');


f_3DPVT = fopen('PVT3D.txt','w');


offset = 1; %121;

% fstart indicates whether Kalman Filter starts
fstart = false;

% fstop indicates whether Kalman Filter stops
fstop = true;

% fwarmup indicates whether Kalman Filter is warming up
fwarmup = false;

% fstart indicates whether Kalman Filter is holding on
fhold = false;

% fcount1 indicates whether Non-NaN count is incremented
fcount1 = false;

% th1 indicates the threshold on which Kalman Filtering should be stopped and reinitialized
th1 = 11;

% count1 notes down how many continuous NaN positions appear
count0 = 0;

% count2 notes down how many continuous non-NaN positions appear
count1 = 0;

th_Ele = 0;

% count_conti counts how many there are continuous measurements in
% total for MHE
count_conti = 0;

for i=offset:N
%     if i==1579
%         fprintf("1422");
%     end
    iValid = find(isfinite(gnssMeas.PrM(i,:)) & gnssMeas.Cn0DbHz(i,:) > 0); %index into valid svid
    svid    = gnssMeas.Svid(iValid)';    
    [gpsEph,iSv] = ClosestGpsEph(allGpsEph,svid,gnssMeas.FctSeconds(i));% full cycle time
    svid = svid(iSv); %svid for which we have ephemeris
    numSvs = length(svid); %number of satellites this epoch
    
    % Satellite elevation filter
    if count1 > 1
        % Compute satellite positions
        prM     = gnssMeas.PrM(i,iValid(iSv))';
        tRx = [ones(numSvs,1)*weekNum(i),gnssMeas.tRxSeconds(i,iValid(iSv))'];
        [svXyzTrx, ~, ~, ~, ~] = ComputeSvXyz(tRx, prM, gpsPvt.allXyzMMMEKF(i-1,:), gpsEph, gpsPvt.allBcMetersEKF(i-1), iono);
        %% Calculate the elevation and azimuth angle
        %Calculate line of sight vectors from user to the satellite
        v = zeros(length(gnssMeas.Svid),3)+NaN;
        v(iValid(iSv),:) = svXyzTrx(:,1:3) - gpsPvt.allXyzMMMEKF(i-1,:);
        
        %Calculate the geodetic latitude and longitude of the user
        llaDegDegM = Xyz2Lla(gpsPvt.allXyzMMMEKF(i-1,:));
        
        %Calculate the rotation matrix to convert an ECEF vector to
        % North, East, Down coordinates, and vice-versa
        RE2N = RotEcef2Ned(llaDegDegM(1),llaDegDegM(2));
        
        %Calculate line of sight vectors from user to the satellite in Ned
        %coordinate system
        vNed = RE2N*v';
        
        % Calculate Elevation angle and azimuth angle
        Ele = asin(-vNed(3,:)./sqrt( sum(vNed.^2) )); % rad
        iValid = find(isfinite(gnssMeas.PrM(i,:)) & gnssMeas.Cn0DbHz(i,:) > 0 & Ele > th_Ele/180*pi);
        gpsPvt.Ele(i,:) = Ele;
        svid    = gnssMeas.Svid(iValid)';    
        [gpsEph,iSv] = ClosestGpsEph(allGpsEph,svid,gnssMeas.FctSeconds(i));% full cycle time
        svid = svid(iSv); %svid for which we have ephemeris
        numSvs = length(svid); %number of satellites this epoch
    end
    
    gpsPvt.iValid(i) = {iValid}; 
    gpsPvt.numSvs(i) = numSvs;
    
    % Initialize count for EKF
    if numSvs >= 4
        % Non-NaN position: continuous locations
        count0 = 0;
        count1 = count1 +1;
        fcount1 = true;
        if Flag_discon(i) == true
            % Continuous locations but discontinuous measurements
            count1 = 1;
            fcount1 = false;
        end
    else % numSvs < 4
        % NaN position: discontinuous locations
        count0 = count0 + 1;
        count1 = 0;
        fcount1 = false;
    end
    
    % Initialize count for MHE
    if Flag_discon(i) == true
        % The current epoch is discontinuous 
        count_conti = 1;
%         % Re-initialize MHE state
%         xo_MHE = [-1509706.868, 6195107.314,149731.63,0,0,0,0,0]';
    else
        % The current epoch is continuous
        count_conti = count_conti + 1;
    end
    
    if count1 > 0
        %% WLS PVT -----------------------------------------------------------------
        %for those svIds with valid ephemeris, pack prs matrix for WlsNav
        prM     = gnssMeas.PrM(i,iValid(iSv))';
%         index_GT = GT_data(:,2) == i;
%         if i > 884 && i < 1036
%         prM = prM - GT_data(index_GT,end-3);
%         end
        
        prSigmaM= gnssMeas.PrSigmaM(i,iValid(iSv))';
        
        prrMps  = gnssMeas.PrrMps(i,iValid(iSv))';
        prrSigmaMps = gnssMeas.PrrSigmaMps(i,iValid(iSv))';
        
        tRx = [ones(numSvs,1)*weekNum(i),gnssMeas.tRxSeconds(i,iValid(iSv))'];
        
        prs = [tRx, svid, prM, prSigmaM, prrMps, prrSigmaMps];
        
        xo(5:8) = zeros(4,1); %initialize speed to zero
        
        
        %% WLS
        [xHat,~,svPos,H,Wpr,Wrr,R_Ac,I_iono_logger,I_trop_logger] = WlsPvt(prs,gpsEph,xo,iono);%compute WLS solution
        xo = xo + xHat;
        
        gpsPvt.SvXyzMMM(iValid(iSv),:,i) = svPos(:, 2:4);
        
        % Calculate the elevation and azimuth angle
        %Calculate line of sight vectors from user to the satellite
        v = zeros(length(gnssMeas.Svid),3)+NaN;
        v(iValid(iSv),:) = svPos(:,2:4) - xo(1:3)';
        
        %Calculate the geodetic latitude and longitude of the user
        llaDegDegM = Xyz2Lla(xo(1:3)');
        
        %Calculate the rotation matrix to convert an ECEF vector to
        % North, East, Down coordinates, and vice-versa
        RE2N = RotEcef2Ned(llaDegDegM(1),llaDegDegM(2));
        
        %Calculate line of sight vectors from user to the satellite in Ned
        %coordinate system
        vNed = RE2N*v';
        
        % Calculate Elevation angle and azimuth angle
        Ele = asin(-vNed(3,:)./sqrt( sum(vNed.^2) )); % rad
        Azi = atan2(vNed(2,:), vNed(1,:)); % rad
       
        % Extract CN0 from gnssMeas
        CN0Sv = gnssMeas.Cn0DbHz(i,iValid(iSv));
        
        % State update of WLS
        % the pseudorange corrected by satellite clock biases
        R = R_Ac + GpsConstants.LIGHTSPEED*svPos(:,5); % R is the pseduorange excluding satellite clock bias, atmospherical delays
        
        gpsPvt.dtsv(i,iValid(iSv)) = svPos(:,5)';
        gpsPvt.dtsvDot(i,iValid(iSv)) = svPos(:,6)';
        gpsPvt.RrMps(i,iValid(iSv)) = svPos(:,7)';
        
        % Pack ECEF positioning results of WLS
        gpsPvt.allXyzMMM(i,:) = xo(1:3);
        
        % Pack lla positioning results of WLS
        gpsPvt.allLlaDegDegM(i,:) = llaDegDegM;
        
        % Pack timing results of WLS
        gpsPvt.allBcMeters(i) = xo(4);
        
        %extract velocity states
        RE2N = RotEcef2Ned(llaDegDegM(1),llaDegDegM(2));
        %NOTE: in real-time code compute RE2N once until position changes
        vNed = RE2N*xo(5:7); %velocity in NED
        gpsPvt.allVelMps(i,:) = vNed;
        gpsPvt.allVelMpsXyz(i,:) = xo(5:7);%velocity in xyz
        gpsPvt.allBcDotMps(i) = xo(8);
        
        gpsPvt.elevation(i,:) = Ele;
        gpsPvt.azimuth(i,:) = Azi;
        
        % extract weights in WLS
        gpsPvt.Wpr(i) = {Wpr};
        
        %% Extended Kalman Filter
        tic
        if (count1 == 1 && fstop == true) || (Flag_discon(i) == true)
            % Warm up Kalman Filter with WLS positioning solutions
            fwarmup = true;
            fstop = false;
            fstart = false;
            fhold = false;
            Xp = [xo(1),xo(5),xo(2),xo(6),xo(3),xo(7),xo(4),xo(8)];
            Pp = zeros(8,8);
            Xhat = [xo(1),xo(5),xo(2),xo(6),xo(3),xo(7),xo(4),xo(8)];
            Phat = zeros(8,8);
        elseif count1 == 2 && fwarmup == true
            % Warm up Kalman Filter with WLS positioning solutions
            fwarmup = true;
            fstop = false;
            fstart = false;
            fhold = false;
            Xhat = gpsPvt.allXhat(i-1,:);
            Phat = gpsPvt.allPhat(:,:,i-1);
            [Xp, Pp] = PredictionEKF(Xhat, Phat, gnssMeas, gpsPvt, i, count1 == 2 && fwarmup == true);
            [Xhat,Phat] = AdjustmentEKF(Xp, Pp, prs,gpsEph,iono);
        elseif (count1 > 2 && fwarmup == true)||(count1 == 1 && fhold == true)||(fstart == true && fcount1 == true)
            % Start Kalman Filter
            fwarmup = false;
            fstop = false;
            fstart = true;
            fhold = false;
            
            Xhat = gpsPvt.allXhat(i-1,:);
            Phat = gpsPvt.allPhat(:,:,i-1);
            [Xp, Pp] = PredictionEKF(Xhat, Phat, gnssMeas, gpsPvt, i, count1 == 2 && fwarmup == true);
            [Xhat,Phat] = AdjustmentEKF(Xp, Pp, prs,gpsEph,iono);
        end
        % Update states of EKF
        % Pack gpsPvt
        gpsPvt.allXp(i,:) = Xp;
        gpsPvt.allPp(:,:,i) = Pp;
        
        gpsPvt.allXhat(i,:) = Xhat;
        gpsPvt.allPhat(:,:,i) = Phat;
        
        % Convert EKF positioning results from ECEF to lla
        llaDegDegMEKF = Xyz2Lla([Xhat(1),Xhat(3),Xhat(5)]);
        
        % Pack lla positioning results of EKF        
        gpsPvt.allLlaDegDegMEKF(i,:) = llaDegDegMEKF;
        
        % Pack ECEF positioning results of EKF
        gpsPvt.allXyzMMMEKF(i,:) = [Xhat(1),Xhat(3),Xhat(5)];
        
        % Compute Pseudorange based on EKF positioning results
        PrMEkf = sqrt(sum((gpsPvt.allXyzMMMEKF(i,:) - svPos(:,2:4))'.^2))';
        
        % Pack timing results of EKF
        gpsPvt.allBcMetersEKF(i) = Xhat(7);
        gpsPvt.allfcMpsEKF(i) = Xhat(8);
        
        % Pack velocity results of EKF     
        gpsPvt.allVelMpsEKF(i,:) = [Xhat(2),Xhat(4),Xhat(6)];%velocity in xyz
        elapsed_time_ekf = elapsed_time_ekf+toc;
        
        
        %% General Update
        A = [i*ones(length(R),1),gnssMeas.FctSeconds(i)*ones(length(R),1),svPos(:,1:5),Ele(iValid(iSv))',Azi(iValid(iSv))',CN0Sv',R,I_iono_logger,I_trop_logger,prM,prSigmaM,prrMps,PrMEkf];
        % Index of time, time, PRN, SvX, SvY, SvZ, dtsv, Elevation, Azimuth, CN0, Pseudorange after corrections, Iono_delay, Trop_delay, Raw Pseudoranges
        fprintf(f_3DPVT,'%d \t %f \t %d \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f\n',A');
        gpsPvt.SvPosR = [gpsPvt.SvPosR; A];
        
        %compute HDOP
        H = [H(:,1:3)*RE2N', ones(numSvs,1)]; %observation matrix in NED
        P = inv(H'*H);%unweighted covariance
        gpsPvt.hdop(i) = sqrt(P(1,1)+P(2,2));
        
        %compute variance of llaDegDegM
        %inside LsPvt the weights are used like this:
        %  z = Hx, premultiply by W: Wz = WHx, and solve for x:
        %  x = pinv(Wpr*H)*Wpr*zPr;
        %  the point of the weights is to make sigma(Wz) = 1
        %  therefore, the variances of x come from  diag(inv(H'Wpr'WprH))
        P = inv(H'*(Wpr'*Wpr)*H); %weighted covariance
        gpsPvt.sigmaLLaM(i,:) = sqrt(diag(P(1:3,1:3)));
        
        %similarly, compute variance of velocity
        P = inv(H'*(Wrr'*Wrr)*H); %weighted covariance
        gpsPvt.sigmaVelMps(i,:) = sqrt(diag(P(1:3,1:3)));
        %%end WLS PVT --------------------------------------------------------------
   elseif (count0 >= th1 && fhold == true)||(count1 == 0 && fstop == true)||(count1 == 0 && fwarmup == true) || (Flag_discon(i) == true)
        % Stop Kalman Filter
        fstop = true;
        fhold = false;
        fstart = false;
        fwarmup = false;
        A = [i,gnssMeas.FctSeconds(i),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
        gpsPvt.SvPosR = [gpsPvt.SvPosR; A];
        gpsPvt.Wpr(i) = {0};
        continue;
    elseif (count1 == 0 && fstart == true)||(count0 < th1 && fhold == true)
        % The current state is NaN
        % but the consecutive NaN measurements are less than the threshold.
        % Hold on Kalman Filter
        tic
        fhold = true;
        fstart = false;
        fwarmup = false;
        fstop = false;
        
        % Just Kalman Prediction
        Xhat = gpsPvt.allXhat(i-1,:);
        Phat = gpsPvt.allPhat(:,:,i-1);
        
        [Xhat, Phat] = PredictionEKF(Xhat, Phat, gnssMeas, gpsPvt, i, count1 == 2 && fwarmup == true);
        
        % Update Kalman Filter
        gpsPvt.allXp(i,:) = Xhat;
        gpsPvt.allPp(:,:,i) = Phat;
        
        gpsPvt.allXhat(i,:) = Xhat;
        gpsPvt.allPhat(:,:,i) = Phat;
        
        gpsPvt.allXyzMMMEKF(i,:) = [Xhat(1),Xhat(3),Xhat(5)];
        gpsPvt.allXyzMMM(i,:) = [Xhat(1),Xhat(3),Xhat(5)];
        
        % Convert EKF positioning results from ECEF to lla
        llaDegDegMEKF = Xyz2Lla([Xhat(1),Xhat(3),Xhat(5)]);
        gpsPvt.allLlaDegDegMEKF(i,:) = llaDegDegMEKF;
        gpsPvt.allLlaDegDegM(i,:) = llaDegDegMEKF;
        
        gpsPvt.allBcMetersEKF(i) = Xhat(7);
        gpsPvt.allBcMeters(i) = Xhat(7);
        
        gpsPvt.allfcMpsEKF(i) = Xhat(8);
        gpsPvt.allBcDotMps(i) = Xhat(8);
        
        gpsPvt.allVelMpsEKF(i,:) = [Xhat(2),Xhat(4),Xhat(6)];
        gpsPvt.allVelMpsXyz(i,:) = gpsPvt.allVelMpsEKF(i,:);
        
        % Compute elapsed time
        elapsed_time_ekf = elapsed_time_ekf+toc;
        
        % Compute satellite positions
        prM     = gnssMeas.PrM(i,iValid(iSv))';
%         index_GT = GT_data(:,2) == i;
%         if i > 884 && i < 1036
%         prM = prM - GT_data(index_GT,end-3);
%         end
        if ~isempty(prM)
            tRx = [ones(numSvs,1)*weekNum(i),gnssMeas.tRxSeconds(i,iValid(iSv))'];
            [svXyzTrx, dtsv, R_Ac, I_iono_logger, I_trop_logger] = ComputeSvXyz(tRx, prM, gpsPvt.allXyzMMMEKF(i,:), gpsEph, gpsPvt.allBcMetersEKF(i), iono);
            
            gpsPvt.SvXyzMMM(iValid(iSv),:,i) = svXyzTrx(:,1:3);
            
            %% Calculate the elevation and azimuth angle
            %Calculate line of sight vectors from user to the satellite
            v = zeros(length(gnssMeas.Svid),3)+NaN;
            v(iValid(iSv),:) = svXyzTrx(:,1:3) - gpsPvt.allXyzMMMEKF(i,:);
            
            %Calculate the geodetic latitude and longitude of the user
            llaDegDegM = Xyz2Lla(gpsPvt.allXyzMMMEKF(i,:));
            
            %Calculate the rotation matrix to convert an ECEF vector to
            % North, East, Down coordinates, and vice-versa
            RE2N = RotEcef2Ned(llaDegDegM(1),llaDegDegM(2));
            
            %Calculate line of sight vectors from user to the satellite in Ned
            %coordinate system
            vNed = RE2N*v';
            
            % Calculate Elevation angle and azimuth angle
            Ele = asin(-vNed(3,:)./sqrt( sum(vNed.^2) )); % rad
            Azi = atan2(vNed(2,:), vNed(1,:)); % rad
            
            %% Extract CN0 from gnssMeas
            CN0Sv = gnssMeas.Cn0DbHz(i,iValid(iSv));
            
            %% State update
            % the pseudorange corrected by satellite clock biases
            R = R_Ac + GpsConstants.LIGHTSPEED*dtsv; % R is the pseduorange excluding satellite clock bias, atmospherical delays
            
            % Compute Pseudorange based on EKF positioning results
            PrMEkf = sqrt(sum((gpsPvt.allXyzMMMEKF(i,:) - svXyzTrx(:,1:3))'.^2))';
            
            prSigmaM= gnssMeas.PrSigmaM(i,iValid(iSv))';
            
            prrMps  = gnssMeas.PrrMps(i,iValid(iSv))';
            
            A = [i*ones(length(R),1),gnssMeas.FctSeconds(i)*ones(length(R),1),svid, svXyzTrx(:,1:3), dtsv, Ele(iValid(iSv))',Azi(iValid(iSv))',CN0Sv',R,I_iono_logger,I_trop_logger,prM,prSigmaM,prrMps,PrMEkf];
            % Index of time, time, PRN, SvX, SvY, SvZ, dtsv, Elevation, Azimuth, CN0, Pseudorange after corrections, Iono_delay, Trop_delay, Raw Pseudoranges
            fprintf(f_3DPVT,'%d \t %f \t %d \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f \t %f\n',A');
            gpsPvt.SvPosR = [gpsPvt.SvPosR; A];           
            gpsPvt.Wpr(i) = {diag(1./prSigmaM)};
        else
            A = [i,gnssMeas.FctSeconds(i),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
            gpsPvt.SvPosR = [gpsPvt.SvPosR; A];
            gpsPvt.Wpr(i) = {0};
        end
    end
    
    %% Moving horizon estimator
    if count_conti == 1 && count1 == 0
        % Clock discontinuity and Satellite discontinuity
        continue;
    else      
        % The window size of MHE, the number of historical data excluding the current epoch
        WindowSize = 8;  
        % N =2 in urban areas
%         WindowSize = 2;
        [xHat_MHE,Sum_numSvs,svid_unique] = MHEstimator(i,gnssMeas,allGpsEph,WindowSize,xo_MHE,weekNum,count_conti,iono);

        if Sum_numSvs<4 || length(svid_unique)<4
            continue;%skip to next epoch
        end

        xo_MHE = xo_MHE + xHat_MHE;
        %     xo_MHE_k_N = xo_MHE_k_N + xHat_k_N;
        %     xo_MHE = xo_MHE_k_N;

        % Convert MHE positioning results from ECEF to lla
        llaDegDegMMHE = Xyz2Lla(xo_MHE(1:3)');

        % Pack ECEF positioning results of MHE
        gpsPvt.allXyzMMMMHE(i,:) = xo_MHE(1:3);

        % Pack lla positioning results of MHE        
        gpsPvt.allLlaDegDegMMHE(i,:) = llaDegDegMMHE;

        % Pack timing results of MHE
        gpsPvt.allBcMetersMHE(i) = xo_MHE(4);
        gpsPvt.allBcDotMpsMHE(i) = xo_MHE(8);       

        %extract velocity states of MHE
        RE2NMHE = RotEcef2Ned(llaDegDegMMHE(1),llaDegDegMMHE(2));
        %NOTE: in real-time code compute RE2N once until position changes
        vNedMHE = RE2NMHE*xo_MHE(5:7); %velocity in NED
        gpsPvt.allVelMpsMHE(i,:) = vNedMHE;
    end
    
end
fclose(f_3DPVT);
time_per_sample_EKF = elapsed_time_ekf/(N-offset+1);
fprintf('\n Processing time per sample of EKF: %s seconds\n', time_per_sample_EKF);

%% RTS Smoothing
% count1 notes down how many continuous non-NaN measurements appear
count1 = 0;
% Backward Kalman Smoothing
tic
for i = N:-1:offset 
    iValid = cell2mat(gpsPvt.iValid(i));
%     iValid = find(isfinite(gnssMeas.PrM(i,:)) & gnssMeas.Cn0DbHz(i,:) > 0); %index into valid svid
    svid    = gnssMeas.Svid(iValid)';
    [~,iSv] = ClosestGpsEph(allGpsEph,svid,gnssMeas.FctSeconds(i));% full cycle time
    
    Xhat = gpsPvt.allXhat(i, :);
    Phat = gpsPvt.allPhat(:,:,i);
    if ~isnan(Xhat(1))
        count1 = count1 +1;
    else
        count1 = 0;            
    end
    
    if i == N
        flag_dis = false;
    else
        flag_dis = Flag_discon(i+1);
    end
    
    if count1 == 0
        % Stop Kalman smoother
        PrMEkfS = 0;
    elseif (count1 == 1) || (flag_dis == true)
        % Warm up Kalman smoother
        % Smoothing starts from the last-time-step state estimated by
        % forward KF
        gpsPvt.allXhatS(i, :) = gpsPvt.allXhat(i, :);
        gpsPvt.allPhatS(:,:,i) = gpsPvt.allPhat(:,:,i);
        
        % Convert EKF positioning results from ECEF to lla
        llaDegDegMEKFS = Xyz2Lla([gpsPvt.allXhatS(i,1),gpsPvt.allXhatS(i,3),gpsPvt.allXhatS(i,5)]);
        
        % Pack lla positioning results of EKF        
        gpsPvt.allLlaDegDegMEKFS(i,:) = llaDegDegMEKFS;
        
        % Pack ECEF positioning results of EKF
        gpsPvt.allXyzMMMEKFS(i,:) = [gpsPvt.allXhatS(i,1),gpsPvt.allXhatS(i,3),gpsPvt.allXhatS(i,5)];
                
        % Pack timing results of EKF
        gpsPvt.allBcMetersEKFS(i) = gpsPvt.allXhatS(i, 7);
        gpsPvt.allfcMpsEKFS(i) = gpsPvt.allXhatS(i, 8);
        
        % Pack velocity results of EKF     
        gpsPvt.allVelMpsEKFS(i,:) = [gpsPvt.allXhatS(i,2),gpsPvt.allXhatS(i,4),gpsPvt.allXhatS(i,6)];%velocity in xyz 
        
        % Compute Pseudorange based on EKF positioning results
        svPos = gpsPvt.SvXyzMMM(iValid(iSv),:,i);
        PrMEkfS = sqrt(sum((gpsPvt.allXyzMMMEKFS(i,:) - svPos())'.^2))';
        
    elseif count1 >= 2 
        % Start Kalman Smoother
        Xs = gpsPvt.allXhatS(i+1, :);
        Ps = gpsPvt.allPhatS(:,:,i+1);
        Xp = gpsPvt.allXp(i+1, :);
        Pp = gpsPvt.allPp(:,:,i+1);
        [Xs, Ps] = SmoothingKF(Xp', Pp, Xhat', Phat, Xs', Ps,i,gnssMeas);
        gpsPvt.allXhatS(i, :) = Xs;
        gpsPvt.allPhatS(:,:,i) = Ps;
        
        % Convert EKF positioning results from ECEF to lla
        llaDegDegMEKFS = Xyz2Lla([Xs(1),Xs(3),Xs(5)]);
        
        % Pack lla positioning results of EKF        
        gpsPvt.allLlaDegDegMEKFS(i,:) = llaDegDegMEKFS;
        
        % Pack ECEF positioning results of EKF
        gpsPvt.allXyzMMMEKFS(i,:) = [Xs(1),Xs(3),Xs(5)];
        
        % Compute Pseudorange based on EKF positioning results
        if isempty(iValid)
            PrMEkfS = 0;
        else
            svPos = gpsPvt.SvXyzMMM(iValid(iSv),:,i);
            PrMEkfS = sqrt(sum((gpsPvt.allXyzMMMEKFS(i,:) - svPos)'.^2))';
        end
        
        % Pack timing results of EKF
        gpsPvt.allBcMetersEKFS(i) = Xs(7);
        gpsPvt.allfcMpsEKFS(i) = Xs(8);
        
        % Pack velocity results of EKF     
        gpsPvt.allVelMpsEKFS(i,:) = [Xs(2),Xs(4),Xs(6)];%velocity in xyz         
    end
    gpsPvt.SvPosRS = [gpsPvt.SvPosRS; PrMEkfS(end:-1:1)];
end
elapsed_time_rts = elapsed_time_ekf+toc;
time_per_sample_RTS = elapsed_time_rts/(N-offset+1);
fprintf('\n Processing time per sample of RTS smoother: %s seconds\n', time_per_sample_RTS);

end %end of function GpsWlsPvt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright 2022 NTUsg.
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

