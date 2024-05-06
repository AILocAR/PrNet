function [bdsEph,iono] = ReadRinex3Nav(fileName)
% [bdsEph,iono] = ReadRinex3Nav(fileName)
%
% Read BDS ephemeris and iono data from an ASCII formatted RINEX 3.03 Nav file.
% Input:
%    fileName - string containing name of RINEX formatted navigation data file
% Output:
% bdsEph: vector of ephemeris data, each element is an ephemeris structure 
%        structure order and orbit variable names follow RINEX 3.0.3 Table A14
% bdsEph(i).PRN     % SV PRN number
% bdsEph(i).Toc     % Time of clock (seconds) BDT
% bdsEph(i).af0     % SV clock bias (seconds)
% bdsEph(i).af1     % SV clock drift (sec/sec)
% bdsEph(i).af2     % SV clock drift rate (sec/sec2)
% bdsEph(i).AODE    % Age of data, ephemeris 
% bdsEph(i).Crs     % Sine harmonic correction to orbit radius (meters)
% bdsEph(i).Delta_n % Mean motion difference from computed value (radians/sec)
% bdsEph(i).M0      % Mean anomaly at reference time (radians)
% bdsEph(i).Cuc     % Cosine harmonic correction to argument of lat (radians)
% bdsEph(i).e       % Eccentricity (dimensionless)
% bdsEph(i).Cus     % Sine harmonic correction to argument of latitude (radians)
% bdsEph(i).Asqrt   % Square root of semi-major axis (meters^1/2)
% bdsEph(i).Toe     % Reference time of ephemeris (seconds)
% bdsEph(i).Cic     % Cosine harmonic correction to angle of inclination (radians)
% bdsEph(i).OMEGA   % Longitude of ascending node at weekly epoch (radians)=OMEGA0 
% bdsEph(i).Cis     % Sine harmonic correction to angle of inclination (radians)
% bdsEph(i).i0      % Inclination angle at reference time (radians)
% bdsEph(i).Crc	    % Cosine harmonic correction to the orbit radius (meters)
% bdsEph(i).omega	% Argument of perigee (radians)
% bdsEph(i).OMEGA_DOT% Rate of right ascension (radians/sec)
% bdsEph(i).IDOT	% Rate of inclination angle (radians/sec)
% bdsEph(i).codeL2  % codes on L2 channel 
% bdsEph(i).BDS_Week % GPS week (to go with Toe), (NOT Mod 1024)
% bdsEph(i).L2Pdata % L2 P data flag
% bdsEph(i).accuracy % SV user range accuracy (meters)
% bdsEph(i).health  % Satellite health
% bdsEph(i).TGD1     % Group delay B1/B3(seconds)
% bdsEph(i).TGD2    % Group delay B2/B3(seconds)
% bdsEph(i).ttx	    % Transmission time of message (seconds)
% bdsEph(i).AODC    %Age of Data Clock
%
% iono: ionospheric parameter structure
%   iono.alpha = [alpha0, alpha1, alpha2, alpha3]
%	iono.beta =  [ beta0,  beta1,  beta2,  beta3]
% if iono data is not present in the Rinex file, iono is returned empty.

fidEph = fopen(fileName);

% Count the number of ephemerides and the lines of header
[numEph,numHdrLines] = countEph(fidEph);

%Now read from the begining again, looking for iono parameters
frewind(fidEph);
iono = readIono(fidEph,numHdrLines);

%initialize ephemeris structure array:
bdsEph = InitializeBdsEph;
bdsEph = repmat(bdsEph,1,numEph);

%now read each ephemeris into bdsEph(j)
%RINEX defines the format in terms of numbers of characters, so that's how we
%read it, e.g. "bdsEph(j).PRN   = str2double(line(1:2));" and so on
for j = 1:numEph
   line         = fgetl(fidEph);
   %The first digit of BDS prn is C,so we start from the second digit
   bdsEph(j).PRN   = str2double(line(2:3));
   %NOTE: we use str2double, not str2double, since str2double handles 'D' for exponent

   %% get Toc (Rinex 3.0 gives this as BDT time yy,mm,dd,hh,mm,ss)
   year   = str2double(line(4:8));%RINEX 3.0 has a 4-digit year

   month  = str2double(line(9:11));
   day    = str2double(line(12:14));
   hour   = str2double(line(15:17));
   minute = str2double(line(18:20));
   second = str2double(line(21:23));
   %convert Toc to gpsTime
   bdsTime      = Utc2Bds([year,month,day,hour,minute,second]);
   bdsEph(j).Toc   = bdsTime(2); %seconds within a week 
   %% get all other parameters
   bdsEph(j).af0   = str2double(line(24:42));%SV clock bias (sec)
   bdsEph(j).af1   = str2double(line(43:61));% SV clock drift (sec/sec)
   bdsEph(j).af2   = str2double(line(62:80));% SV clock drift rate (sec/sec2)
   
   line = fgetl(fidEph);
   bdsEph(j).AODE  = str2double(line(5:23));
   bdsEph(j).Crs   = str2double(line(24:42));
   bdsEph(j).Delta_n = str2double(line(43:61));
   bdsEph(j).M0    = str2double(line(62:80));
   
   line = fgetl(fidEph);
   bdsEph(j).Cuc   = str2double(line(5:23));
   bdsEph(j).e     = str2double(line(24:42));
   bdsEph(j).Cus   = str2double(line(43:61));
   bdsEph(j).Asqrt = str2double(line(62:80));

   line=fgetl(fidEph);
   bdsEph(j).Toe   = str2double(line(5:23));
   bdsEph(j).Cic   = str2double(line(24:42));
   bdsEph(j).OMEGA = str2double(line(43:61));
   bdsEph(j).Cis   = str2double(line(62:80));

   line = fgetl(fidEph); 
   bdsEph(j).i0        =  str2double(line(5:23));
   bdsEph(j).Crc       = str2double(line(24:42));
   bdsEph(j).omega     = str2double(line(43:61));
   bdsEph(j).OMEGA_DOT = str2double(line(62:80));
   
   line = fgetl(fidEph);
   bdsEph(j).IDOT      = str2double(line(5:23));
%    bdsEph(j).codeL2    = str2double(line(24:42));
   bdsEph(j).BDS_Week  = str2double(line(43:61));
%    bdsEph(j).L2Pdata   = str2double(line(62:80));

   line = fgetl(fidEph);
   bdsEph(j).accuracy  = str2double(line(5:23));
   bdsEph(j).health    = str2double(line(24:42));
   bdsEph(j).TGD1       = str2double(line(43:61));
   bdsEph(j).TGD2      = str2double(line(62:80));
   
   line = fgetl(fidEph);
   bdsEph(j).ttx           = str2double(line(5:23));
   bdsEph(j).AODC  = str2double(line(24:42));
end
fclose(fidEph);

end %end of function ReadRinexNav
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Count ephemerides
function [numEph,numHdrLines] = countEph(fidEph,fileName) 
%utility function for ReadRinexNav
%Read past the header, and then read to the end, counting ephemerides:
numHdrLines = 0;
bFoundHeader = false;

%Find "END OF HEADER"
while ~bFoundHeader  %Read past the header
    numHdrLines = numHdrLines+1;
    line = fgetl(fidEph);
    if ~ischar(line), break, end
    k = strfind(line,'END OF HEADER');
    if ~isempty(k),
        bFoundHeader = true;
        break
    end
end
if ~bFoundHeader
    error('Error reading file: %s\nExpected RINEX header not found',fileName);
end
%Now HEADER ends. Start to read to the end of the file
numEph = -1;% because the end of the file, a NaN line, is also counted
while 1
    numEph = numEph+1;
    line = fgetl(fidEph);
    if line == -1 %the end of the file
        break;  
%     elseif length(line)~=80 %Each line in a BDS ephemeris file has 80 columns
%         %use this opportunity to check line is the right length
%         %because in the rest of ReadRinexNav we depend on line being this length
%         error('Incorrect line length encountered in RINEX file'); 
    end
end
%check that we read the expected number of lines: 8 lines per ephemeris
if rem(numEph,8)~=0
  error('Number of nav lines in %s should be divisible by 8',fileName);
end
numEph = numEph/8; %8 lines per ephemeris

end %end of function countEph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function iono = readIono(fidEph,numHdrLines)
%utility function to read thru the header lines, and find iono parameters

iono = []; %return empty if iono not found
bIonoAlpha=false; 
% bIonoBeta=false;

for i = 1:numHdrLines, 
    line = fgetl(fidEph); 
        % Look for iono parameters, and read them in
    if ~isempty(strfind(line,'IONOSPHERIC CORR')) && ~isempty(strfind(line,'BDSA')) %line contains iono alpha parameters
        ii = strfind(line,'IONOSPHERIC CORR');
        iono.alpha=str2double(line(5:ii-7));
        %If we have 4 parameters then we have the complete iono alpha
        bIonoAlpha = (length(iono.alpha)==4);
    end
    if ~isempty(strfind(line,'IONOSPHERIC CORR')) && ~isempty(strfind(line,'BDSB'))%line contains the iono beta parameters
        ii = strfind(line,'IONOSPHERIC CORR');
        iono.beta=str2double(line(5:ii-7));
        %If we have 4 parameters then we have the complete iono beta
        bIonoBeta = (length(iono.beta)==4);
    end
end

if ~(bIonoAlpha && bIonoBeta)
   %if we didn't get both alpha and beta iono correctly, then return empty iono
   iono=[];
end

end %end of function readIono
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function bdsEph = InitializeBdsEph
%utility function to initialize the ephemeris structure
bdsEph.PRN         = 0;
bdsEph.Toc         = 0;
bdsEph.af0         = 0;
bdsEph.af1         = 0;
bdsEph.af2         = 0;
bdsEph.AODE        = 0;
bdsEph.Crs         = 0;
bdsEph.Delta_n     = 0;
bdsEph.M0          = 0;
bdsEph.Cuc         = 0;
bdsEph.e           = 0;
bdsEph.Cus         = 0;
bdsEph.Asqrt       = 0;
bdsEph.Toe         = 0;
bdsEph.Cic         = 0;
bdsEph.OMEGA       = 0;
bdsEph.Cis         = 0;
bdsEph.i0          = 0;
bdsEph.Crc         = 0;
bdsEph.omega       = 0;
bdsEph.OMEGA_DOT   = 0;
bdsEph.IDOT        = 0;
% bdsEph.codeL2      = 0;%bdsEph.spare
bdsEph.BDS_Week    = 0;
% bdsEph.L2Pdata     = 0;%bdsEph.spare
bdsEph.accuracy    = 0;
bdsEph.health      = 0;
bdsEph.TGD1         = 0;
bdsEph.TGD2        = 0;%
bdsEph.ttx         = 0;%
bdsEph.AODC= 0;%

end %end of function InitializeEph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Copyright 2016 Google Inc.
% Copyright 2021 Nanyang Technological University
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
