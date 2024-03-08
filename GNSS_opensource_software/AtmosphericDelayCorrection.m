function [I_iono, I_trop] = AtmosphericDelayCorrection(svXyzTrx,xyz0,iono,tRx)

%Calculate line of sight vectors from user to the satellite
v = svXyzTrx - xyz0;

%Calculate the geodetic latitude and longitude of the user
llaDegDegM = Xyz2Lla(xyz0);

%Calculate the rotation matrix to convert an ECEF vector to 
% North, East, Down coordinates, and vice-versa
 RE2N = RotEcef2Ned(llaDegDegM(1),llaDegDegM(2));

%Calculate line of sight vectors from user to the satellite in Ned
%coordinate system
vNed = RE2N*v';

%% 0. Calculate Elevation angle and azimuth angle 
E = asin(-vNed(3)/sqrt( sum(vNed.^2) )); % rad
A = atan2(vNed(2), vNed(1)); % rad

%% 1.Calculate the earth's central angle ? between the user position and the earth projection of ionospheric pierce (intersection) point
Phi = 0.0137/(E/pi+0.11) - 0.022; % semi-circles

%% 2.Calculate the geodetic latitude ?i of the earth projection of the ionospheric pierce (intersection) point (semi-circles).
phi_i = llaDegDegM(1)/180 + Phi*cos(A); % semi-circles
if phi_i > 0.416
    phi_i = 0.416;
elseif phi_i < -0.416
    phi_i = -0.416; 
end

%% 3.Calculate the geodetic longitude of the earth projection of the ionospheric pierce (intersection) point (semi-circles).
lamda_i = llaDegDegM(2)/180 + Phi*sin(A)/cos(phi_i*pi);% semi-circles

%% 4.Calculate the geomagnetic latitude of the earth projection of the ionospheric pierce (intersection) point (mean ionospheric height assumed 350 km) (semi-circles)
phi_m = phi_i + 0.064*cos((lamda_i-1.617)*pi); % semi-circles

%% 5.Calculate amplitude (AMP) and period (PER) of Klobuchar Model with the ionospheric parameters ? and ?. 
AMP = iono.alpha(1) + iono.alpha(2)*phi_m + iono.alpha(3)*phi_m^2 + iono.alpha(4)*phi_m^3;
if AMP < 0
    AMP = 0;
end

PER = iono.beta(1) + iono.beta(2)*phi_m + iono.beta(3)*phi_m^2 + iono.beta(4)*phi_m^3;
if PER < 72000
    PER = 72000;
end

%% 6.Calculate the local time
t = 43200*lamda_i + tRx(2) - 86400 * floor(tRx(2)/86400); % sec

if t >= 86400
    t = t-86400;
elseif t < 0
    t = t+86400;
end

%% 7.Calculate the phase x of Klobuchar Model
x = 2*pi*(t-50400)/PER; % radians

%% 8.Calculate the obliquity factor F (dimensionless)
F = 1 + 16*(0.53-E/pi)^3;

%% 9.Calculate the ionospheric delay
if abs(x) < 1.57
    I = F * (5e-9 + AMP*(1 - x^2/2 + x^4/24));
else
    I = F * 5e-9;
end
I_iono = I * GpsConstants.LIGHTSPEED;

%% Calculate tropospheric delay
if llaDegDegM(3) > 200000.0
    mFactor = 0.0;
elseif llaDegDegM(3) < 0
    mFactor = 1.0;
else
    mFactor=exp(-1*llaDegDegM(3)*1.33e-4);
end
    
I_trop = mFactor * 2.47/(sin(E)+0.0121); 
end