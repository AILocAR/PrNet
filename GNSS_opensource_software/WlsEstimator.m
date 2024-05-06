function [xo, Ele, Azi, CN0Sv, R, llaDegDegM] = WlsEstimator(xo, gnssMeas, gpsEph, iValid, iSv, weekNum, i, svid, numSvs, iono)

%for those svIds with valid ephemeris, pack prs matrix for WlsNav
prM     = gnssMeas.PrM(i,iValid(iSv))';
prSigmaM= gnssMeas.PrSigmaM(i,iValid(iSv))';
prrMps  = gnssMeas.PrrMps(i,iValid(iSv))';
prrSigmaMps = gnssMeas.PrrSigmaMps(i,iValid(iSv))';

tRx = [ones(numSvs,1)*weekNum(i),gnssMeas.tRxSeconds(i,iValid(iSv))'];

prs = [tRx, svid, prM, prSigmaM, prrMps, prrSigmaMps];

xo(5:8) = zeros(4,1); %initialize speed to zero


%% WLS
[xHat,~,svPos,H,Wpr,Wrr,R_Ac,I_iono_logger,I_trop_logger] = WlsPvt(prs,gpsEph,xo,iono);%compute WLS solution

xo = xo + xHat;

%% Calculate the elevation and azimuth angle
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

%% Extract CN0 from gnssMeas
CN0Sv = gnssMeas.Cn0DbHz(i,iValid(iSv));

%% State update
% the pseudorange corrected by satellite clock biases
R = R_Ac + GpsConstants.LIGHTSPEED*svPos(:,5); % R is the pseduorange excluding satellite clock bias, atmospherical delays












end