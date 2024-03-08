function Rz = RotZ(theta)

% Rotation matrix for rotation around Z axis 
% North, East, Down coordinates, and vice-versa
%
% inputs: theta (degrees)
% output: Re2n,   3x3 unitary rotation matrix = 
%              [cos(theta),  sin(theta),  0;
%               -sin(theta), cos(theta),  0;
%               0,           0,           1]
%
% Example:   vR = RotZ*v 


%Author: Xu Weng
%Open Source code for PrNet and Outdoor AR

%CHECK INPUTS
if any(size(theta)~=[1,1])
    error('Inputs theta must be scalars')
end

D2R = pi/180; %degrees to radians scale factor
thetaRad=D2R*theta(:); 

ctheta = cos(thetaRad);
stheta = sin(thetaRad);

Rz = zeros(3,3);
Rz(1,1) = ctheta;
Rz(1,2) = stheta;
Rz(1,3) = 0;

Rz(2,1) = -stheta;
Rz(2,2) = ctheta;
Rz(2,3) = 0;

Rz(3,1) = 0;
Rz(3,2) = 0;
Rz(3,3) = 1;

end %end of function RotZ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%