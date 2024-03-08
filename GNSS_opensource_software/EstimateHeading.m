function Hat_Heading = EstimateHeading(gpsPvt)
% Estimate the heading of the smartphone
% By Xu Weng @NTU,Sg

% The number of time steps
N = length(gpsPvt.FctSeconds);


Hat_Heading = zeros(N, 3);

for i= 1: N
    % Calculate the unit vector from the current epoch to the next epoch
    % Convert the position at the current epoch from ECEF to LLA
    if i == N
        v_xyz = gpsPvt.allXyzMMMEKF(N,:) - gpsPvt.allXyzMMMEKF(N-1,:);    
    else
        v_xyz = gpsPvt.allXyzMMMEKF(i+1,:) - gpsPvt.allXyzMMMEKF(i,:);      
    end
    
    % Normalized v_xyz
    v_xyz = v_xyz/norm(v_xyz);
    
    % Convert v_xyz from ECEF to NED
    % 1. Convert the position at the current epoch from ECEF to LLA
    current_lla = Xyz2Lla(gpsPvt.allXyzMMMEKF(i,:));
    
    % 2. Calculate the rotation matrix to convert an ECEF vector to
    %    North, East, Down coordinates, and vice-versa   
    RE2N = RotEcef2Ned(current_lla(1), current_lla(2));
    
    % 3. Convert v_xyz to v_ned
    v_ned = RE2N*v_xyz';
    
    Hat_Heading(i, :) = v_ned';    
end

end