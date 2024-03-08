function [Xhat,Phat] = AdjustmentEKF(Xp, Pp, prs,gnssEph,iono)
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

%
% Outputs: xHat: estimate of state update


%Author: Frank van Diggelen
%Modified by Xu Weng
%Open Source code for processing Android GNSS Measurements

jWk=1; jSec=2; jSv=3; jPr=4; jPrSig=5; jPrr=6; jPrrSig=7;%index of columns

[bOk,numVal] = checkInputs(prs, gnssEph, Xp);
if ~bOk
    error('inputs not right size, or not properly aligned with each other')
end

Xhat=[]; Phat=[]; 
xyz0 = [Xp(1), Xp(3), Xp(5)];
bc = Xp(7);
fc = Xp(8);

if numVal<4
  return
end
ttxWeek = prs(:,jWk); %week of tx. Note - we could get a rollover, when ttx_sv
%goes negative, and it is handled in GpsEph2Pvt, where we work with fct
ttxSeconds =  prs(:,jSec) - prs(:,jPr)/GpsConstants.LIGHTSPEED; %ttx by sv clock 
% this is accurate satellite time of tx, because we use actual pseudo-ranges 
% here, not corrected ranges
% write the equation for pseudorange to see the rx clock error exactly cancel
% to get precise GPS time: we subtract the satellite clock error from sv time, 
% as done next:
dtsv = GpsEph2Dtsv(gnssEph,ttxSeconds);
dtsv = dtsv(:); %make into a column for compatibility with other time vectors
ttx = ttxSeconds - dtsv; %subtract dtsv from sv time to get true gps time

%calculate satellite position at ttx
[svXyzTtx,dtsv,svXyzDot,dtsvDot]=GpsEph2Pvt(gnssEph,[ttxWeek,ttx]);

svXyzTrx = svXyzTtx; %initialize svXyz at time of reception

% Eliminate the effect on the orbit of satellites of the rotation of the earch  
for i=1:length(gnssEph)
        % calculate tflight from, bc and dtsv
        dtflight = (prs(i,jPr)-bc)/GpsConstants.LIGHTSPEED + dtsv(i);
        % Use of bc: bc>0 <=> pr too big <=> tflight too big.
        %   i.e. trx = trxu - bc/GpsConstants.LIGHTSPEED
        % Use of dtsv: dtsv>0 <=> pr too small <=> tflight too small.
        %   i.e ttx = ttxsv - dtsv
        svXyzTrx(i,:) = FlightTimeCorrection(svXyzTtx(i,:), dtflight);
        % Ionospheric Delay Correction
        [I_iono,I_trop] = AtmosphericDelayCorrection(svXyzTrx(i,:),xyz0,iono,[prs(i,jWk), prs(i,jSec)]); % meter
        prs(i,jPr) = prs(i,jPr) - I_iono - I_trop;
end

%calculate line of sight vectors and ranges from satellite to xo
v = xyz0(:)*ones(1,numVal,1) - svXyzTrx';%v(:,i) = vector from sv(i) to xyz0
%Normalization
range = sqrt( sum(v.^2) );
v = v./(ones(3,1)*range); % line of sight unit vectors from sv to xo

% Initialize the measurement matrix
C = zeros(2*numVal,8);

% Initialize the pseudorange rate measurement matrix
rrMps = zeros(numVal,1);

% Initialize the measurement residuals  matrix
b = zeros(2*numVal,1);

% Initialize the measurement noise covariance matrix
Er = zeros(2*numVal,1);

for i = 1:numVal
    
    % Construct the measurement matrix
    C(2*i-1:2*i,:)=[v(1,i), 0,      v(2,i), 0,      v(3,i), 0,      1, 0;
                    0,      v(1,i), 0     , v(2,i), 0,      v(3,i), 0, 1];   
                          
    % Construct the measurement residuals
    % Pseduorange measurement
    b(2*i-1,1) = prs(i,jPr) - (range(i) + bc - GpsConstants.LIGHTSPEED*dtsv(i));
    
    % Pseduorange rate measurement
    %range rate = [satellite velocity] dot product [los from sv to xo]
    rrMps(i) = svXyzDot(i,:)*v(:,i);% v(n)*l(n)
    Vp =[Xp(2),Xp(4),Xp(6)];
    b(2*i,1) = prs(i,jPrr) + GpsConstants.LIGHTSPEED*dtsvDot(i) + rrMps(i) - (Vp*v(:,i) + fc);
    
    % Construct the measurement noise covariance matrix
    Er(2*i-1,1) = prs(i,jPrSig)^2;
    Er(2*i,1) = prs(i,jPrrSig)^2;
end

R = diag(Er);

% Calculate Kalman Gain
K = Pp*C'*(C*Pp*C'+R)^(-1);

% Calculate the posterior estimation
Xhat = Xp + K*b;

% Update estimation covariance
Phat = (eye(8)-K*C)*Pp;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [bOk,numVal] = checkInputs(prs, gpsEph, Xp)
%utility function for WlsPvt
jWk=1; jSec=2; jSv=3; jPr=4; jPrSig=5; jPrr=6; jPrrSig=7;%index of columns

bOk=false;
%check inputs
numVal=size(prs,1);  
if (max(prs(:,jSec))-min(prs(:,jSec)))> eps
  return
elseif length(gpsEph)~=numVal
    return
elseif any(prs(:,jSv) ~= [gpsEph.PRN]')
    return
elseif  any(size(Xp) ~= [8,1])
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
% Copyright 2022 Nanyang Technological University
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
