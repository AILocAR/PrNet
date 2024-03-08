function [xHat_k,z,svPos,H,Wpr_rr] = WlsPvtMHE2ndOrder(prs_stack,gpsEph_stack,xo,iono,time_i)
% [xHat,z,svPos,H,Wpr,Wrr] = WlsPvt(prs,gpsEph,xo)
% calculate a weighted least squares PVT solution, xHat
% given pseudoranges, pr rates, and initial state
%
% Inputs:
%  prs: matrix of raw pseudoranges, and pr rates, each row of the form:

%  [trxWeek,trxSeconds,sv,prMeters,prSigmaMeters,prrMps,prrSigmaMps]
%   trxWeek, trxSeconds: Rx time of measurement
%      where trxSeconds = seconds in the current GPS week

%   sv: satellite id number
%   prMeters, prSigmaMeters: pseudorange and standard deviation (meters)
%   prrMps, prrSigmaMps: pseudorange rate and standard deviation (m/s)
%   gpsEph: matching vector of GPS ephemeris struct, defined in ReadRinexNav

%   xo: initial (previous) state, [x,y,z,bc,xDot,yDot,xDot,bcDot]'
%       in ECEF coordinates(meters and m/s)
%
% Outputs: xHat: estimate of state update
%          z = [zPr; zPrr] a-posteriori residuals (measured-calculated)

%          svPos: matrix of calculated sv positions and sv clock error:
%                 [sv prn, x,y,z (ecef m), dtsv (s),xDot,yDot,zDot, dtsvDot]

%          H: H observation matrix corresponding to svs in svPos(:,1)

%          Wpr,Wrr Weights used in WlsPvt = 1/diag(sigma measurements)
%                  use these matrices to compute variances of xHat
%
% The PVT solution = xo + xHat, in ECEF coordinates

% For unweighted solution, set all sigmas = 1

%Author: Modified by Xu Weng and Minyi Lin based on Frank van Diggelen's codes
%Open Source code for processing Android GNSS Measurements

%% Initialization
xHat=[]; z=[]; H=[]; svPos=[];

% Initialize the states at k-N epoch by xo
xyz0 = xo(1:3);
bc = xo(4);
V0 = xo(5:7);
fc = xo(8);
X_k = [xyz0(1),V0(1),xyz0(2),V0(2),xyz0(3),V0(3),bc,fc]'; % State at the kth epoch


jWk=1; jSec=2; jSv=3; jPr=4; jPrSig=5; jPrr=6; jPrrSig=7;%index of columns

%prs =[tRx, svid, prM, prSigmaM, prrMps, prrSigmaMps];

% Window Size
L = length(prs_stack);

% Allocate space to initial state each epoch in the window
HatX_k = zeros(8,L); 
% HatX_k = zeros(7,L); 


%% Calculate the satellite positions for each epoch in the moving horizon window

% Allocate space to satellite positions and dtsv for each epoch in the moving horizon window
svXyzTrx_cell = cell([L,1]);
svXyzTtx_cell = cell([L,1]);
svXyzDot_cell = cell([L,1]);
dtsv_cell = cell([L,1]);
dtsvDot_cell = cell([L,1]);
v_cell = cell([L,1]);
diag_Wpr_rr = [];


% Total number of available ephemeris values in the current window
Sum_numVal = 0;

% Log number of available sv each epoch in the current window
numVal_MHE  = zeros(L,1);

for k=L:-1:1
    
    prs = cell2mat(prs_stack(k));
    gpsEph = cell2mat(gpsEph_stack(k));
    [bOk,numVal] = checkInputs(prs, gpsEph, xo);
    numVal_MHE(k) = numVal;
    Sum_numVal = Sum_numVal+numVal;
    
    if ~bOk
        error('inputs not right size, or not properly aligned with each other')
    end
        
    %jWk=1; jSec=2; jSv=3; jPr=4; jPrSig=5; jPrr=6; jPrrSig=7
    ttxWeek = prs(:,jWk); %week of tx. Note - we could get a rollover, when ttx_sv
    %goes negative, and it is handled in GpsEph2Pvt, where we work with fct
        
    ttxSeconds =  prs(:,jSec) - prs(:,jPr)/GpsConstants.LIGHTSPEED; %ttx by sv clock
    
    % this is accurate satellite time of tx, because we use actual pseudo-ranges
    % here, not corrected ranges
    % write the equation for pseudorange to see the rx clock error exactly cancel
    % to get precise GPS time: we subtract the satellite clock error from sv time,
    % as done next:
    dtsv = GpsEph2Dtsv(gpsEph,ttxSeconds);
    dtsv = dtsv(:); %make into a column for compatibility with other time vectors
    ttx = ttxSeconds - dtsv; %subtract dtsv from sv time to get true gps time
    
    %calculate satellite position at ttx
    [svXyzTtx,dtsv,svXyzDot,dtsvDot] = GpsEph2Pvt(gpsEph,[ttxWeek,ttx]);
    
    svXyzTrx = svXyzTtx; %initialize svXyz at time of reception
    svXyzTtx_cell(k) = {svXyzTtx};
    svXyzTrx_cell(k) = {svXyzTrx};
    dtsv_cell(k) = {dtsv};
    dtsvDot_cell(k) = {dtsvDot};
    svXyzDot_cell(k) = {svXyzDot};
    
    %Compute weights ---------------------------------------------------
    for i = 1:numVal
       diag_Wpr_rr = [diag_Wpr_rr;prs(i,jPrSig);prs(i,jPrrSig)]; 
    end
    
end

Wpr_rr = diag(1./diag_Wpr_rr);


if Sum_numVal<4
    return
end

%% WLS
%iterate on this next part till change in pos & line of sight vectors converge
xHat=zeros(8,1);
% xHat=zeros(7,1);
dx=xHat+inf;
whileCount=0; maxWhileCount=1000;
%we expect the while loop to converge in < 10 iterations, even with initial
%position on other side of the Earth (see Stanford course AA272C "Intro to GPS")
while norm(dx) > GnssThresholds.MAXDELPOSFORNAVM
    whileCount=whileCount+1;
    assert(whileCount < maxWhileCount,...
        'while loop did not converge after %d iterations',whileCount);
    % Allocate space to Packed measurements
    bb = [];

    % Allocate space to Packed Geometry Matrix
    CC = [];
    
    % Initialize state update matrix
    AA = eye(8);
    for j=L:-1:1
        
        % Fetch the jth measurements from package cells
        prs = cell2mat(prs_stack(j));
        dtsv = cell2mat(dtsv_cell(j));
        svXyzTtx = cell2mat(svXyzTtx_cell(j));
        svXyzTrx = cell2mat(svXyzTrx_cell(j));
        svXyzDot = cell2mat(svXyzDot_cell(j));
        dtsvDot = cell2mat(dtsvDot_cell(j));
        numVal = numVal_MHE(j);
        
        % Calculate the initial state each epoch
        if j==L
            HatX_k(:,j) = X_k;
        else
            % Define the state transition matrix
            % Calculate the sampling interval first
            prs_l = cell2mat(prs_stack(j+1));
            Ts = (prs_l(1,jSec)+prs_l(1,jWk)*GpsConstants.WEEKSEC)-(prs(1,jSec)+prs(1,jWk)*GpsConstants.WEEKSEC);
            
            a0 = zeros(2,2);
            a = [1, Ts; 0, 1];
            a1 = [1, 0; 0, 1];
            A = [a, a0, a0, a0;
                a0, a, a0, a0;
                a0, a0, a, a0;
                a0, a0, a0, a];
%             A = [a, a0, a0, a0;
%                 a0, a, a0, a0;
%                 a0, a0, a, a0;
%                 0,  0,  0, 1];
            AA = A\AA;
            HatX_k(:,j) = A\HatX_k(:,j+1);
        end
        
        % Assign estimated states at the jth epoch
        xyz_k = HatX_k(1:2:5,j);
        bc_k = HatX_k(7,j);     
        V_k = HatX_k(2:2:6,j);
        fc_k = HatX_k(8,j);
                
        for i=1:numVal
            % calculate tflight from, bc and dtsv
            dtflight = (prs(i,jPr)-bc_k)/GpsConstants.LIGHTSPEED + dtsv(i);
            % Use of bc: bc>0 <=> pr too big <=> tflight too big.
            %   i.e. trx = trxu - bc/GpsConstants.LIGHTSPEED
            % Use of dtsv: dtsv>0 <=> pr too small <=> tflight too small.
            %   i.e ttx = ttxsv - dtsv
            svXyzTrx(i,:) = FlightTimeCorrection(svXyzTtx(i,:), dtflight);
            % Ionospheric Delay Correction
            [I_iono,I_trop] = AtmosphericDelayCorrection(svXyzTrx(i,:),xyz_k',iono,[prs(i,jWk), prs(i,jSec)]); % meter
            prs(i,jPr) = prs(i,jPr) - I_iono - I_trop;
        end
        
        %calculate line of sight vectors and ranges from satellite to xo
        v = xyz_k*ones(1,numVal,1) - svXyzTrx'; %v(:,i) = vector from sv(i) to xyz0 
        
        range = sqrt( sum(v.^2) ); %Normolize v
        v = v./(ones(3,1)*range); % line of sight unit vectors from sv to xo ...corresponding to ak(n)
        
        % Initialize the measurement matrix
        C = zeros(2*numVal,8);
        
        % Initialize the pseudorange rate measurement matrix
        rrMps = zeros(numVal,1);
        
        % Initialize the measurement residuals  matrix
        b = zeros(2*numVal,1);
        
        for i = 1:numVal
            
            % Construct the measurement matrix
            C(2*i-1:2*i,:)=[v(1,i), 0,      v(2,i), 0,      v(3,i), 0,      1, 0;
                            0,      v(1,i), 0     , v(2,i), 0,      v(3,i), 0, 1];
            
            % Construct the measurement residuals
            % Pseduorange measurement
            b(2*i-1,1) = prs(i,jPr) - (range(i) + bc_k - GpsConstants.LIGHTSPEED*dtsv(i));
            
            % Pseduorange rate measurement
            %range rate = [satellite velocity] dot product [los from sv to xo]
            rrMps(i) = svXyzDot(i,:)*v(:,i);% v(n)*a(n)

            b(2*i,1) = prs(i,jPrr) + GpsConstants.LIGHTSPEED*dtsvDot(i) + rrMps(i) - (V_k'*v(:,i) + fc_k);
            
        end
        
        
        % Store v in the corresponding cell
        v_cell(j) = {v};
        
        svPos=[prs(:,3),svXyzTrx,dtsv(:)];
        
        %% packaging b and C and then calculate estimate
        bb = [bb;b];
        
        CC = [CC;C*AA];
    end
    %z = Hx, premultiply by W: Wz = WHx, and solve for x:
    % dx = pinv(Wpr*CC)*Wpr*bb;
    dx = pinv(Wpr_rr*CC)*Wpr_rr*bb;
        
    % State Update 
    xHat = xHat+dx;% delta[x,vx,y,vy,z,vz,bc,fc]'
    xyz0 = xyz0(:) + dx(1:2:5);
    bc = bc + dx(7);
    V0 = V0 + dx(2:2:6);
    fc = fc + dx(8);
    
    X_k = [xyz0(1),V0(1),xyz0(2),V0(2),xyz0(3),V0(3),bc,fc]'; % Estimated State at the kth epoch
end
% the state update at the kth epoch
xHat_k = [xHat(1),xHat(3),xHat(5),xHat(7),xHat(2),xHat(4),xHat(6),xHat(8)]'; 

% % the state update at the k-Nth epoch
% xHat = AA*xHat;
% xHat = [xHat(1),xHat(3),xHat(5),xHat(7),xHat(2),xHat(4),xHat(6),xHat(8)]'; %[x,y,z,bc,vx,vy,vz,fc]

end %end of function WlsPvt


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [bOk,numVal] = checkInputs(prs, gpsEph, xo)
%utility function for WlsPvt
jWk=1; jSec=2; jSv=3; jPr=4; jPrSig=5; jPrr=6; jPrrSig=7;%index of columns

bOk=false;
%check inputs
numVal=size(prs,1);

if (max(prs(:,jSec))-min(prs(:,jSec)))> eps
    return    
elseif length(gpsEph)~=numVal
    return
elseif any(prs(:,jSv) ~= [gpsEph.PRN]')  %PRN（pseudorandom noise)
    %an important element of code division multiple access (CDMA) based satellite navigation systems.
    %Each satellite within a GNSS constellation has a UNIQUE PRN code that it transmits as part of the C/A navigation message.
    %This code allows any receiver to identify exactly which satellite(s) it is receiving.
    return
elseif  any(size(xo) ~= [8,1])
    return
elseif size(prs,2)~=7
    return
else
    bOk = true;
end

%We insist that gpsEph and prs are aligned first.
%ClosestGpsEph.m does this, and passes back indices for prs - this is the way to
%do it right, so we don't have nested searches for svId

end %end of function checkInputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright 2016 Google Inc.
% Copyright 2022 NTU, Singapore
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



