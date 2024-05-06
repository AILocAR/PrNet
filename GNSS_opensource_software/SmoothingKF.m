function [Xs, Ps] = SmoothingKF(Xp, Pp, Xhat, Phat, Xs, Ps, iepoch, gnssMeas)
% State X = [prM, prrMps]
% Input:
%      Xp --- Prior estimation at k+1 time step
%      Pp --- Prior covariance matrix at k+1 time step
%      Xhat --- Posterior estimation at k time step
%      Phat --- Posterior covariance matrix at k time step
%      Xs --- Smoothed estimation at k+1 time step
%      Ps --- Posterior covariance matrix of the estimation at k+1 time step

% Output:
%      Xs --- Smoothed estimation at k time step
%      Ps --- Posterior covariance matrix of the estimation at k time step
% Author: Xu Weng @ NTUsg

% For low speed movement, Full Time in Seconds is precise enough


Ts = gnssMeas.FctSeconds(iepoch+1) - gnssMeas.FctSeconds(iepoch);

% State transition matrix 
a0 = zeros(2,2);
a = [1, Ts; 0, 1];
A = [a, a0, a0, a0;
     a0, a, a0, a0;
     a0, a0, a, a0;
     a0, a0, a0, a];

% Smoothing Gain
G = Phat*A'*Pp^(-1);

Xs = Xhat+G*(Xs - Xp);

Ps = Phat + G*(Ps - Pp)*G';

end

% Copyright 2022 Nanyang Technological University
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
% http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.