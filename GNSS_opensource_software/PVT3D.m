clc;
clear;

%load the measurements including time, satellite PRN, satellite positions
%and the meausured pseudorange corrected by satellite clock biases

file_3D_PVT = load('PVT3D.txt');
t_index = 1; % the epoch index
f_row_index = 1;% the row index of file_3D_PVT
[row_f,col_f] = size(file_3D_PVT);

xHat = zeros(4,1); %initial state of the approximate position of the receiver: [center of the Earth, tu=0]'

while f_row_index <= row_f
    i = 1;% the first sv at the t_index epoch
    while file_3D_PVT(f_row_index,1) == t_index && f_row_index <= row_f
       SvPVT(i,:) = file_3D_PVT(f_row_index,3:7);
       i = i+1;  
       f_row_index = f_row_index+1;
    end
    
    %% calculate the receiver's position at the t_index epoch
    
    %calculate line of sight vectors from the approximate position to the
    %satellites
    [numVal,~] = size(SvPVT);
    xyz0 = xHat(1:3);%initialize xyz
    tu = xHat(4);%initialize receiver clock bias
    DeltaX = xHat+inf;
    DeltaX0 = xHat+0;
    count = 0;
    %Least Squres method
    while norm(DeltaX) > 20
        count = count+1;
        %calculate unit vectors from the user to satellites
        a = SvPVT(:,2:4)'- xyz0(:)*ones(1,numVal);%a(:,i) = vector from xyz0 to sv(i)
        %Normalization
        RHat = sqrt( sum(a.^2) );
        a = a./(ones(3,1)*RHat); % line of sight unit vectors from xyz0 to sv
        H = [a', ones(numVal,1)]; % H matrix = [unit vector,1]
        
        %calculate the range residual
        DeltaR = RHat(:) - (SvPVT(:,5) - tu);
        
        %calculate the position offset
        DeltaX0 = DeltaX;
        DeltaX = pinv(H)*DeltaR;
        
        % update xyz0 and tu
        %xHat=xHat+dx;
        xyz0=xyz0(:)+DeltaX(1:3);
        tu = tu-DeltaX(4);
    end
    UserPVT(t_index,:) = [xyz0',tu];
    xHat = [xyz0',tu];
    t_index = t_index+1;
end