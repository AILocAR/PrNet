function csvNmeaFileName = ReadNmeaFile(dirName, fileName)
% Read NMEA standard file
% Return GPS time and positioning results in a .csv file

% Make extended file name
if dirName(end)~='/'
    dirName = [dirName,'/']; %add /
end

extendedFileName = [dirName,fileName];
fprintf('\nReading file %s\n',extendedFileName)

% Create the .csv file
csvNmeaFileName = [dirName,fileName(1:end-4), 'csv'];
csvfileID = fopen(csvNmeaFileName,'w');

% Open the NMEA file
txtfileID = fopen(extendedFileName,'r');
% Whether the file is opened successfully
if txtfileID<0
    error('file ''%s'' not found',extendedFileName);
end

% Read the NMEA file
line = fgetl(txtfileID);

while ischar(line)
    
    if contains(line,'$GPGGA,')      
        %Now 'line' contains the GPGGA message
        line = strrep(line,'$GPGGA,',''); %remove '$GPGGA'; 
        line = strrep(line,'N,',''); %remove 'N'; 
        line = strrep(line,'W,',''); %remove 'W'; 
        line = strrep(line,'M,',''); %remove 'M'; 
        line = strrep(line,'*',''); %remove '*';
        line = regexprep(line,'[a-zA-Z]',''); %remove English characteristics;
        fprintf(csvfileID,'%s\n',line);
    end
    line = fgetl(txtfileID);

end
fclose(txtfileID);
fclose(csvfileID);
end