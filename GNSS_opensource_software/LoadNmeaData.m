function GroundTruth = LoadNmeaData(dirName, fileName)
% Extract the date of the data file
indexlist = strfind(dirName, '/');
index = indexlist(end-2);
% index = indexlist(end);
Year = str2double(dirName(index+1:index+4));
Month = str2double(dirName(index+6:index+7));
Day = str2double(dirName(index+9:index+10));

% Read ground truth file
fprintf('\nReading file %s\n',fileName)

% Load raw ground truth data
GroundTruth0 = load(fileName);
[GT_row, ~] = size(GroundTruth0);

% Create space for output
GroundTruth = zeros(GT_row, 9);

for i = 1:GT_row
    
    % Extract time
    Hour = floor(GroundTruth0(i,1)/10000);
    Minute = floor((GroundTruth0(i,1)-Hour*10000)/100);
    Second = GroundTruth0(i,1)-Hour*10000 - Minute*100;
    
    % Check a jump of a day
    % The groundtruth is 1Hz        
    if i ~= 1 && Hour == 0 && Minute == 0 && floor(Second) ==0 && ~contains(fileName, "10Hz")
        Day = Day+1;
    end
    % The groundtruth is 10Hz
    if i ~= 1 && Hour == 0 && Minute == 0 && Second ==0 && contains(fileName, "10Hz")
        Day = Day+1;
    end
          
    % What about a jump of month and year???
    
    utcTime = [Year,Month,Day,Hour,Minute,Second];
    [gpsTime, ~] = Utc2Gps(utcTime);
    
    % Extract latitude, longitude and amplitude
    GT_latitude = floor(GroundTruth0(i,2)/100) + (GroundTruth0(i,2) - floor(GroundTruth0(i,2)/100)*100)/60;
    GT_longitude = floor(GroundTruth0(i,3)/100) + (GroundTruth0(i,3) - floor(GroundTruth0(i,3)/100)*100)/60;
    GT_Amplitude = GroundTruth0(i,7) + GroundTruth0(i,8); % Height = height relative to the mean sea level + geoidal seperation
    
    GroundTruthXyz = Lla2Xyz([GT_latitude, -GT_longitude, GT_Amplitude]);
    GroundTruth(i,:) = [i, gpsTime,-GT_longitude, GT_latitude, GT_Amplitude, GroundTruthXyz];
    
end

end