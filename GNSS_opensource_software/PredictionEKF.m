function [Xp, Pp] = PredictionEKF(Xhat, Phat, gnssMeas, gpsPvt, iepoch, flag2)
% State X = [x, vx, y,vy, z, vz, deltat, deltaf]
% Input:
%      Xhat --- Posterior estimation of state of last epoch
%      Pe --- Posterior covariance matrix
%      gnssMeas --- GNSS raw measurement
%      iepoch --- Index of current epoch
%      flag2 --- flag2 indicates if the current epoch is the second epoch
%      in warm-up process
% Output:
%      Xp --- Prior estimation of current state
%      Pp --- Prior covariance matrix

% For low speed movement, Full Time in Seconds is precise enough
% Author: Xu Weng
% if iepoch < 3
%     error('Kalman Filtering should start at the third epoch');
% end
Ts = gnssMeas.FctSeconds(iepoch) - gnssMeas.FctSeconds(iepoch - 1);

% Power spectrum density of clock bias noise
% St = ((gnssMeas.BiasUncertaintyNanos(iepoch)*1e-9*GpsConstants.LIGHTSPEED)^2+(gnssMeas.BiasUncertaintyNanos(iepoch-1)*1e-9*GpsConstants.LIGHTSPEED)^2)/Ts^2;
% Sf = ((gnssMeas.DriftUncertaintyNanosPerSecond(iepoch)*1e-9)^2+(gnssMeas.DriftUncertaintyNanosPerSecond(iepoch-1)*1e-9)^2)/Ts^2;

% Power spectrum density of velocity noise (acceleration) on three axes
% Power spectrum density of clock drift noise
% if iepoch == 3
%     Svxyz = ((gpsPvt.allVelMpsXyz(iepoch-1,:)-gpsPvt.allVelMpsXyz(iepoch-2,:))/Ts).^2;
%     Sf = ((gpsPvt.allBcDotMps(iepoch-1)-gpsPvt.allBcDotMps(iepoch-2))/Ts).^2;
%     St = ((gpsPvt.allBcMeters(iepoch-1)-gpsPvt.allBcMeters(iepoch-2))/Ts - gpsPvt.allBcDotMps(iepoch-1)).^2;
% elseif iepoch == 4
%     Svxyz = ((gpsPvt.allVelMpsEKF(iepoch-1,:)-gpsPvt.allVelMpsXyz(iepoch-2,:))/Ts).^2;
%     Sf = ((gpsPvt.allfcMpsEKF(iepoch-1)-gpsPvt.allBcDotMps(iepoch-2))/Ts).^2;
%     St = ((gpsPvt.allBcMetersEKF(iepoch-1)-gpsPvt.allBcMeters(iepoch-2))/Ts - gpsPvt.allfcMpsEKF(iepoch-1)).^2;
% elseif iepoch > 4
%     Svxyz = ((gpsPvt.allVelMpsEKF(iepoch-1,:)-gpsPvt.allVelMpsEKF(iepoch-2,:))/Ts).^2;
%     Sf = ((gpsPvt.allfcMpsEKF(iepoch-1)-gpsPvt.allfcMpsEKF(iepoch-2))/Ts).^2;
%     St = ((gpsPvt.allBcMetersEKF(iepoch-1)-gpsPvt.allBcMetersEKF(iepoch-2))/Ts - gpsPvt.allfcMpsEKF(iepoch-1)).^2;
% end

if flag2 == true
    Svxyz = ((gpsPvt.allVelMpsXyz(iepoch,:)-gpsPvt.allVelMpsEKF(iepoch-1,:))/Ts).^2;
    Sf = ((gpsPvt.allBcDotMps(iepoch)-gpsPvt.allfcMpsEKF(iepoch-1))/Ts).^2;
    St = ((gpsPvt.allBcMeters(iepoch)-gpsPvt.allBcMetersEKF(iepoch-1))/Ts - gpsPvt.allBcDotMps(iepoch)).^2;
else
    Svxyz = ((gpsPvt.allVelMpsEKF(iepoch-1,:)-gpsPvt.allVelMpsEKF(iepoch-2,:))/Ts).^2;
    Sf = ((gpsPvt.allfcMpsEKF(iepoch-1)-gpsPvt.allfcMpsEKF(iepoch-2))/Ts).^2;
    St = ((gpsPvt.allBcMetersEKF(iepoch-1)-gpsPvt.allBcMetersEKF(iepoch-2))/Ts - gpsPvt.allfcMpsEKF(iepoch-1)).^2;
end

Svx = Svxyz(1);
Svy = Svxyz(2);
Svz = Svxyz(3);

% Covariance matrix of process noise
Q0 = zeros(2,2);
Qx = [Svx*Ts^3/3, Svx*Ts^2/2; Svx*Ts^2/2, Svx*Ts];
Qy = [Svy*Ts^3/3, Svy*Ts^2/2; Svy*Ts^2/2, Svy*Ts];
Qz = [Svz*Ts^3/3, Svz*Ts^2/2; Svz*Ts^2/2, Svz*Ts];
Qt = [St*Ts+Sf*Ts^3/3, Sf*Ts^2/2; Sf*Ts^2/2, Sf*Ts];

Qe = [Qx, Q0, Q0, Q0;
     Q0, Qy, Q0, Q0;
     Q0, Q0, Qz, Q0;
     Q0, Q0, Q0, Qt];

% State transition matrix 
a0 = zeros(2,2);
a = [1, Ts; 0, 1];
A = [a, a0, a0, a0;
     a0, a, a0, a0;
     a0, a0, a, a0;
     a0, a0, a0, a];
 
Xp = A*Xhat';
Pp = A*Phat*A'+Qe;

end

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