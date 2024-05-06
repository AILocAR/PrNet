function [svXyzTrx, dtsv, R_Ac, I_iono_logger, I_trop_logger] = ComputeSvXyz(tRx, PrM, xyz, gnssEph, bc, iono)

% Input:
% tRx      ---- the time when signals are received
% PrM      ---- Pseudorange 
% xyz      ---- user's position at tRx
% gnssEph  ---- GNSS ephemeris 
% bc       ---- user clock bias
% iono     ---- parameters describing ionospheric delays

% Output:
% svXyzTrx ---- Positions of satellites at tRx 
% R_Ac     ---- Pseudorange measurements after atmospheric corrections
% I_iono_logger ----- Ionospheric delays
% I_trop_logger ----- Tropospheric delays


ttxWeek = tRx(:,1); %week of tx. Note - we could get a rollover, when ttx_sv
%goes negative, and it is handled in GpsEph2Pvt, where we work with fct

ttxSeconds =  tRx(:,2) - PrM/GpsConstants.LIGHTSPEED; %ttx by sv clock 
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

R_Ac = zeros(length(gnssEph),1)+NaN;
I_iono_logger = zeros(length(gnssEph),1)+NaN;
I_trop_logger = zeros(length(gnssEph),1)+NaN;

% Ionospheric Delay Correction
for i = 1:length(gnssEph)
    % calculate tflight from, bc and dtsv
    dtflight = (PrM(i)-bc)/GpsConstants.LIGHTSPEED + dtsv(i);
    
    svXyzTrx(i,:) = FlightTimeCorrection(svXyzTtx(i,:), dtflight);
    
    [I_iono,I_trop] = AtmosphericDelayCorrection(svXyzTrx(i,:),xyz,iono,tRx(i,:)); % meter
    PrM(i) = PrM(i) - I_iono - I_trop;
    I_iono_logger(i) = I_iono;
    I_trop_logger(i) = I_trop;
    R_Ac(i) = PrM(i);
end


end