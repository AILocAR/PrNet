function Hat_Wls_dBc = EstimateWlsDbc(gpsPvt, PrmErr)
% Estimate WLS estimation error of user clock bias
% By Xu Weng @NTU,Sg

% The number of time steps
N = length(gpsPvt.FctSeconds);


Hat_Wls_dBc = zeros(N,1);

for i= 1: N
    % Get the satellite at the current time step
    iValid = gpsPvt.SvPosR(:, 1)==i;
    
    % Get the pseudorange measurements excluding satellite clock bias, atmospherical delays 
    R = gpsPvt.SvPosR(iValid, end-3);
        
    % Get the visible satellite positions
    SvXyz = gpsPvt.SvPosR(iValid, 4:6);
    
    % Get the receiver position estimation at the current time step
    WlsXyz = gpsPvt.allXyzMMM(i, :)';
        
    % Calculate Geometry Matrix
    % 0. Number of visible satellites
    numVal = length(R);
    if numVal == 1
        continue;
    end
    % 1. Calculate line of sight vectors from satellite to xo
    v = WlsXyz*ones(1,numVal,1) - SvXyz';%v(:,i) = vector from sv(i) to xyz0
    % 2. Calculate ranges from satellite to xo
    range = sqrt( sum(v.^2) );
    % 3. Normalization
    v = v./(ones(3,1)*range); % line of sight unit vectors from sv to xo
    % 4. Construct geometry matrix    
    H = [v', ones(numVal,1)]; % H matrix = [unit vector,1]
    
    % Calculate the pseudo-inverse of H
    Wpr = cell2mat(gpsPvt.Wpr(i));
    pinvH = pinv(Wpr*H)*Wpr;
       
    % Get the pseudorange measurement errors at the current time step            
    PrErrPerStep = PrmErr(iValid);
    
    % Compute estimation error of user clock bias
    Hat_Wls_dBc(i) = pinvH(end,:)*PrErrPerStep;
    
end

end