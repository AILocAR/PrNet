function [gpsEph,iono] = ReadRinex3NavGps(fileName)
% [gpsEph,iono] = ReadRinexNav(fileName)
%
% Read GPS ephemeris and iono data from an ASCII formatted RINEX 2.10 Nav file.
% Input:
%    fileName - string containing name of RINEX formatted navigation data file
% Output:
% gpsEph: vector of ephemeris data, each element is an ephemeris structure 
%        structure order and orbit variable names follow RINEX 2.1 Table A4
%        clock variable names af0, af1, af2 follow IS GPS 200
% gpsEph(i).PRN     % SV PRN number
% gpsEph(i).Toc     % Time of clock (seconds)
% gpsEph(i).af0     % SV clock bias (seconds)
% gpsEph(i).af1     % SV clock drift (sec/sec)
% gpsEph(i).af2     % SV clock drift rate (sec/sec2)
% gpsEph(i).IODE    % Issue of data, ephemeris 
% gpsEph(i).Crs     % Sine harmonic correction to orbit radius (meters)
% gpsEph(i).Delta_n % Mean motion difference from computed value (radians/sec)
% gpsEph(i).M0      % Mean anomaly at reference time (radians)
% gpsEph(i).Cuc     % Cosine harmonic correction to argument of lat (radians)
% gpsEph(i).e       % Eccentricity (dimensionless)
% gpsEph(i).Cus     % Sine harmonic correction to argument of latitude (radians)
% gpsEph(i).Asqrt   % Square root of semi-major axis (meters^1/2)
% gpsEph(i).Toe     % Reference time of ephemeris (seconds)
% gpsEph(i).Cic     % Cosine harmonic correction to angle of inclination (radians)
% gpsEph(i).OMEGA   % Longitude of ascending node at weekly epoch (radians)
% gpsEph(i).Cis     % Sine harmonic correction to angle of inclination (radians)
% gpsEph(i).i0      % Inclination angle at reference time (radians)
% gpsEph(i).Crc	    % Cosine harmonic correction to the orbit radius (meters)
% gpsEph(i).omega	% Argument of perigee (radians)
% gpsEph(i).OMEGA_DOT% Rate of right ascension (radians/sec)
% gpsEph(i).IDOT	% Rate of inclination angle (radians/sec)
% gpsEph(i).codeL2  % codes on L2 channel 
% gpsEph(i).GPS_Week % GPS week (to go with Toe), (NOT Mod 1024)
% gpsEph(i).L2Pdata % L2 P data flag
% gpsEph(i).accuracy % SV user range accuracy (meters)
% gpsEph(i).health  % Satellite health
% gpsEph(i).TGD     % Group delay (seconds)
% gpsEph(i).IODC    % Issue of Data, Clock 
% gpsEph(i).ttx	    % Transmission time of message (seconds)
% gpsEph(i).Fit_interval %fit interval (hours), zero if not known
%
% iono: ionospheric parameter structure
%   iono.alpha = [alpha0, alpha1, alpha2, alpha3]
%	iono.beta =  [ beta0,  beta1,  beta2,  beta3]
% if iono data is not present in the Rinex file, iono is returned empty.

fidEph = fopen(fileName);
[numEph,numHdrLines] = countEph(fidEph);

%Now read from the begining again, looking for iono parameters
frewind(fidEph);
iono = readIono(fidEph,numHdrLines);

%initialize ephemeris structure array:
gpsEph = InitializeGpsEph;
gpsEph = repmat(gpsEph,1,numEph);

%now read each ephemeris into gpsEph(j)
%RINEX defines the format in terms of numbers of characters, so that's how we
%read it, e.g. "gpsEph(j).PRN   = str2num(line(1:2));" and so on
for j = 1:numEph
   line         = fgetl(fidEph);
   %The first digit of GPS prn is G,so we start from the second digit
   gpsEph(j).PRN   = str2num(line(2:3));
   %NOTE: we use str2num, not str2double, since str2num handles 'D' for exponent

   %% get Toc (Rinex 3.0 gives this as GPS time yy,mm,dd,hh,mm,ss)
   year   = str2double(line(4:8));%RINEX 3.0 has a 4-digit year

   month  = str2double(line(9:11));
   day    = str2double(line(12:14));
   hour   = str2double(line(15:17));
   minute = str2double(line(18:20));
   second = str2double(line(21:23));
   %convert Toc to gpsTime
   gpsTime      = Utc2Gps([year,month,day,hour,minute,second]);
   gpsEph(j).Toc   = gpsTime(2);
   %% get all other parameters
   gpsEph(j).af0   = str2num(line(24:42));
   gpsEph(j).af1   = str2num(line(43:61));
   gpsEph(j).af2   = str2num(line(62:80));
   
   line = fgetl(fidEph);
   gpsEph(j).IODE  = str2num(line(5:23));
   gpsEph(j).Crs   = str2num(line(24:42));
   gpsEph(j).Delta_n = str2num(line(43:61));
   gpsEph(j).M0    = str2num(line(62:80));
   
   line = fgetl(fidEph);
   gpsEph(j).Cuc   = str2num(line(5:23));
   gpsEph(j).e     = str2num(line(24:42));
   gpsEph(j).Cus   = str2num(line(43:61));
   gpsEph(j).Asqrt = str2num(line(62:80));

   line=fgetl(fidEph);
   gpsEph(j).Toe   = str2num(line(5:23));
   gpsEph(j).Cic   = str2num(line(24:42));
   gpsEph(j).OMEGA = str2num(line(43:61));
   gpsEph(j).Cis   = str2num(line(62:80));

   line = fgetl(fidEph); 
   gpsEph(j).i0        =  str2num(line(5:23));
   gpsEph(j).Crc       = str2num(line(24:42));
   gpsEph(j).omega     = str2num(line(43:61));
   gpsEph(j).OMEGA_DOT = str2num(line(62:80));
   
   line = fgetl(fidEph);
   gpsEph(j).IDOT      = str2num(line(5:23));
   gpsEph(j).codeL2    = str2num(line(24:42));
   gpsEph(j).GPS_Week  = str2num(line(43:61));
   gpsEph(j).L2Pdata   = str2num(line(62:80));

   line = fgetl(fidEph);
   gpsEph(j).accuracy  = str2num(line(5:23));
   gpsEph(j).health    = str2num(line(24:42));
   gpsEph(j).TGD       = str2num(line(43:61));
   gpsEph(j).IODC      = str2num(line(62:80));
   
   line = fgetl(fidEph);
   gpsEph(j).ttx           = str2num(line(5:23));
   gpsEph(j).Fit_interval  = str2num(line(24:42));
end
fclose(fidEph);

end %end of function ReadRinexNav
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [numEph,numHdrLines] = countEph(fidEph,fileName)
%utility function for ReadRinexNav
%Read past the header, and then read to the end, counting ephemerides:
numHdrLines = 0;
bFoundHeader = false;
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
%Now read to the end of the file
numEph = -1;
while 1
    numEph = numEph+1;
    line = fgetl(fidEph);
    if line == -1 %the end of the file
        break;  
%     elseif length(line)~=79
%         %use this opportunity to check line is the right length
%         %because in the rest of ReadRinexNav we depend on line being this length
%         error('Incorrect line length encountered in RINEX file'); 
    end
end
%check that we read the expected number of lines:
if rem(numEph,8)~=0
  error('Number of nav lines in %s should be divisible by 8',fileName);
end
numEph = numEph/8; %8 lines per ephemeris

end %end of function countEph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function iono = readIono(fidEph,numHdrLines)
%utility function to read thru the header lines, and find iono parameters

iono = []; %return empty if iono not found
bIonoAlpha=false; bIonoBeta=false;

for i = 1:numHdrLines, 
    line = fgetl(fidEph); 
        % Look for iono parameters, and read them in
    if ~isempty(strfind(line,'ION ALPHA')) %line contains iono alpha parameters
        ii = strfind(line,'ION ALPHA');
        iono.alpha=str2num(line(1:ii-1));
        %If we have 4 parameters then we have the complete iono alpha
        bIonoAlpha = (length(iono.alpha)==4);
    end
    if ~isempty(strfind(line,'ION BETA'))%line contains the iono beta parameters
        ii = strfind(line,'ION BETA');
        iono.beta=str2num(line(1:ii-1));
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


function gpsEph = InitializeGpsEph
%utility function to initialize the ephemeris structure
gpsEph.PRN         = 0;
gpsEph.Toc         = 0;
gpsEph.af0         = 0;
gpsEph.af1         = 0;
gpsEph.af2         = 0;
gpsEph.IODE        = 0;
gpsEph.Crs         = 0;
gpsEph.Delta_n     = 0;
gpsEph.M0          = 0;
gpsEph.Cuc         = 0;
gpsEph.e           = 0;
gpsEph.Cus         = 0;
gpsEph.Asqrt       = 0;
gpsEph.Toe         = 0;
gpsEph.Cic         = 0;
gpsEph.OMEGA       = 0;
gpsEph.Cis         = 0;
gpsEph.i0          = 0;
gpsEph.Crc         = 0;
gpsEph.omega       = 0;
gpsEph.OMEGA_DOT   = 0;
gpsEph.IDOT        = 0;
gpsEph.codeL2      = 0;
gpsEph.GPS_Week    = 0;
gpsEph.L2Pdata     = 0;
gpsEph.accuracy    = 0;
gpsEph.health      = 0;
gpsEph.TGD         = 0;
gpsEph.IODC        = 0;
gpsEph.ttx         = 0;
gpsEph.Fit_interval= 0;

end %end of function InitializeEph
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Copyright 2016 Google Inc.
% 
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
