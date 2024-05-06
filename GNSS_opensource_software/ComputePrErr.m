function [PrErr, Wls_dBc, h_end] = ComputePrErr(gpsPvt, GroundTruth, tRx_i)
% Compute pseudorange errors
% By Xu Weng @NTU,Sg

% The number of time steps
N = length(gpsPvt.FctSeconds);
% PrErr = zeros(length(gpsPvt.SvPosR(:, 1)))+NaN;
PrErr = [];
h_end = [];

Wls_dBc = zeros(N,1);

for i= tRx_i: N
    % Get the satellite at the current time step
    iValid = gpsPvt.SvPosR(:, 1)==i;
    
    % Get the pseudorange measurements excluding satellite clock bias, atmospherical delays 
    R = gpsPvt.SvPosR(iValid, 11);
    
    % Get the receiver clock bias estimation at the current time step
    bc_hat = gpsPvt.allBcMeters(i); 
    
    % Get the visible satellite positions
    SvXyz = gpsPvt.SvPosR(iValid, 4:6);
    
    % Get ground truth of receiver
    if size(GroundTruth, 1) == 1
        GTXyz = GroundTruth';
    else
        GTXyz = GroundTruth(i - tRx_i + 1,:)';
    end
    
    % Calculate Geometry Matrix
    % 0. Number of visible satellites
    numVal = length(R);
    if numVal == 1
        PrErrPerStep = 0;
        PrErr = [PrErr; PrErrPerStep];
        h_end = [h_end;0];
        continue;
    end
    % 1. Calculate line of sight vectors from satellite to xo
    v = GTXyz*ones(1,numVal,1) - SvXyz';%v(:,i) = vector from sv(i) to xyz0
    % 2. Calculate ranges from satellite to xo
    range = sqrt( sum(v.^2) );
    % 3. Normalization
    v = v./(ones(3,1)*range); % line of sight unit vectors from sv to xo
    % 4. Construct geometry matrix    
    H = [v', ones(numVal,1)]; % H matrix = [unit vector,1]
    
    % Calculate the pseudo-inverse of H
    Wpr = cell2mat(gpsPvt.Wpr(i));
    pinvH = pinv(Wpr*H)*Wpr;
    
    % Use the last row of pinvH to calculate the estimation error of
    % receiver clock error
    % Then, use least squares method to calculate the pseudorange errors
    G = eye(numVal)-ones(numVal,1)*pinvH(end,:);
    h_end = [h_end; pinvH(end,:)'];
    % Calculate pseudorange residuals
    % Pseudorange residuals = pseudorange noise error + estimation error of
    % receiver clock bias
    PrRes = R -range'- bc_hat;
        
    PrErrPerStep = pinv(G)*PrRes;
    PrErr = [PrErr; PrErrPerStep];
    
    % Compute estimation error of user clock bias
    Wls_dBc(i) = pinvH(end,:)*PrErrPerStep;
    
end

end