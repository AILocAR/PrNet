clc;
clear;
%ProcessGnssMeasScript.m, script to read GnssLogger output, compute and plot:
% pseudoranges, C/No, and weighted least squares PVT solution
%
% you can run the data in pseudoranges log files provided for you: 

% prFileName = 'Pixel5_GnssLog.txt';
% prFileName = 'SamsungS20Ultra_GnssLog.txt';
prFileName = 'Pixel4_GnssLog.txt';
% prFileName = 'gnss_log.txt';

%Groundtruth file name
% gtNmeaFileName = 'SPAN_Pixel5_10Hz.nmea';
% gtNmeaFileName = 'SPAN_SamsungS20Ultra_10Hz.nmea';
gtNmeaFileName = 'SPAN_Pixel4_10Hz.nmea';
% gtNmeaFileName = 'span_log.nmea';

% as follows
% 1) copy everything from GitHub google/gps-measurement-tools/ to 
%    a local directory on your machine
% 2) change 'dirName = ...' to match the local directory you are using:
dirName ='../Data/RouteU/Testing/GSDC/2021-04-28-US-SJC-1/Pixel4/supplemental';
% 3) run ProcessGnssMeasScript.m script file 
param.llaTrueDegDegM = [];

indexlist = strfind(dirName, '/');

% index1 = indexlist(end);
% filename_date = dirName(index1+1:end);

index1 = indexlist(end-2);
index2 = indexlist(end-1);
filename_date = dirName(index1+1:index2-1);

%Author: Frank van Diggelen
%Modified by Aaron (Xu Weng)
%Open Source code for processing Android GNSS Measurements

%% parameters
%param.llaTrueDegDegM = [];
%enter true WGS84 lla, if you know it:
param.llaTrueDegDegM = [37.422578, -122.081678, -28];%Charleston Park Test Site

%% Set the data filter and Read log file
dataFilter = SetDataFilter; 
% FullBiasNanos                 FullBiasNanos ~= 0; 
% ConstellationType             ConstellationType==1 (GPS);
% State                         bitand(State,2^0) & bitand(State,2^3)

%% Read Groundtruth data
csvNmeaFileName = ReadNmeaFile(dirName, gtNmeaFileName);
GroundTruth_10Hz = LoadNmeaData(dirName, csvNmeaFileName);
% GroundTruth  = GroundTruth_10Hz(1:10:end,:);
GroundTruth  = GroundTruth_10Hz;
figure_gt = figure;
plot(GroundTruth(1:10:end,4),GroundTruth(1:10:end,5),'color','#008744','linewidth',6);hold on;
% plot(GroundTruth(round(length(GroundTruth)/2):end,4),GroundTruth(round(length(GroundTruth)/2):end,5),'color','#008744','linewidth',6);hold on;
% plot(GroundTruth(2286:2314,4),GroundTruth(2286:2314,5),'color','#008744','linewidth',6);hold on;
% legend("Ground Truth");
xlabel("Longitude",'FontSize',40);
ylabel("Latitude",'FontSize',40);
axis tight;
hold on;
% axis([-122.5, -122.0, 37.3, 37.8]);


% SJC area
I = imread([dirName '/SJC.png']);
axis([-121.902, -121.88, 37.325, 37.34]);
h = image('XData',[-121.902, -121.88],'YData',[37.34,37.325],'CData',I);%note the latitude (y-axis) is flipped in vertical direction

% % MTV area
% I = imread([dirName '/MTV.png']);
% axis([-122.5, -121.8, 37.3, 37.7]);
% h = image('XData',[-122.5, -121.8],'YData',[37.7,37.3],'CData',I);%note the latitude (y-axis) is flipped in vertical direction

uistack(h,'bottom'); %move the image to the bottom of current stack
saveas(figure_gt,[dirName,'/Results_GT.fig']);

% [gnssRaw,gnssAnalysis] = ReadGnssLogger(dirName,prFileName,dataFilter);
[gnssRaw,gnssAnalysis] = ReadGnssLoggerNew(dirName,prFileName,dataFilter);
if isempty(gnssRaw), return, end

%% Read Analysis results from Google's Analysis Software
% AnFileName = 'test1_gnss_log_2022_06_24_14_30_16_derived.txt';
% gnssDerived = ReadGnssAnalysis(dirName,AnFileName,dataFilter);


%% Get online ephemeris from Nasa ftp, first compute UTC Time from gnssRaw:
fctSeconds = 1e-3*double(gnssRaw.allRxMillis(end)); %Full cycle time in seconds
utcTime = Gps2Utc([],fctSeconds);
[allGpsEph, iono] = GetNasaHourlyEphemeris(utcTime,dirName);
% allGpsEph = GetCsnoHourlyEphemeris(utcTime,dirName);
if isempty(allGpsEph), return, end

%% process raw measurements, compute pseudoranges:
% [gnssMeas] = ProcessGnssMeas(gnssRaw);
[gnssMeas] = ProcessGnssMeasNew(gnssRaw);

%% plot pseudoranges and pseudorange rates
h1 = figure;
[colors] = PlotPseudoranges(gnssMeas,prFileName);
h2 = figure;
PlotPseudorangeRates(gnssMeas,prFileName,colors);
h3 = figure;
PlotCno(gnssMeas,prFileName,colors);

%% Detect Smartphone discontinuity
% Detect Smartphone clock discontinuity
Threshold_discontinuity = 10;
Flag_discon1 = [0;diff(gnssMeas.FctSeconds)] > Threshold_discontinuity;

% Detect Smartphone PrM discontinuity
Flag_discon2 = sum(abs([zeros(1,size(gnssMeas.PrM, 2));diff(gnssMeas.PrM)]) > 5e4, 2) > 0;

% Reported HW clock discontinuouty
Flag_discon3 = [0;diff(gnssMeas.ClkDCount)~=0];%binary, 1 <=> clock discontinuity

% Combine 2 discontinuity
Flag_discon = Flag_discon1 | Flag_discon2 | Flag_discon3;

%% compute WLS position and velocity
gpsPvt = GpsWlsPvtEKF_test(gnssMeas,allGpsEph,iono,Flag_discon);

f_3D_WLS_PVT = fopen('WLS_PVT_3D.csv','w');
fprintf(f_3D_WLS_PVT,'%f , %f , %f \n',gpsPvt.allLlaDegDegM');
fclose(f_3D_WLS_PVT);

%% Mountain View
figure;
plot(gpsPvt.allLlaDegDegM(:,2),gpsPvt.allLlaDegDegM(:,1),'m','linewidth',2);hold on;
plot(gpsPvt.allLlaDegDegMEKF(:,2),gpsPvt.allLlaDegDegMEKF(:,1),'c','linewidth',2);hold on;
legend("WLS Results","EKF");
axis tight;
hold on;

% SFO area
% I = imread([dirName '/SFO.png']);
% h = image('XData',[-122.35, -122.1],'YData',[37.5,37.4],'CData',I);%note the latitude (y-axis) is flipped in vertical direction % SFO

% SJC area
I = imread([dirName '/SJC.png']); 
axis([-121.902, -121.88, 37.325, 37.34]);
h = image('XData',[-121.902, -121.88],'YData',[37.34,37.325],'CData',I);%note the latitude (y-axis) is flipped in vertical direction

% % MTV area
% I = imread([dirName '/MTV.png']); 
% axis([-122.5, -121.8, 37.3, 37.7]);
% % axis([-122.5, -122.0, 37.3, 37.8]);
% h = image('XData',[-122.5, -121.8],'YData',[37.7,37.3],'CData',I);%note the latitude (y-axis) is flipped in vertical direction
% % h = image('XData',[-122.5, -122.0],'YData',[37.7,37.4],'CData',I);%note the latitude (y-axis) is flipped in vertical direction

uistack(h,'bottom'); %move the image to the bottom of current stack

%% Comparison with ground truth receiver
[GT_row, GT_clo] = size(GroundTruth);
[tRx_row, tRx_clo] = size(gnssMeas.tRxSeconds);
tRx_i = 1;
tRx_j = 1;

% The groundtruth is later than Android raw measurements
if GroundTruth(1,3) - gnssMeas.tRxSeconds(tRx_i, find(~isnan(gnssMeas.tRxSeconds(tRx_i, :)), 1)) > 0
    for i = tRx_i:tRx_row
        if abs(GroundTruth(1,3)-gnssMeas.tRxSeconds(i, find(~isnan(gnssMeas.tRxSeconds(i, :)), 1))) <= 0.5
            % Align the first epoch
            break;
        end
    end
    tRx_i = i;
    % The ground truth file is 10 Hz
    GroundTruth_cache = GroundTruth(1:10:GT_row,:);  
%     % The ground truth file is 1 Hz
%     GroundTruth_cache = GroundTruth(1:GT_row,:);
    
% The groundtruth is earlier than Android raw measurements
elseif GroundTruth(1,3) - gnssMeas.tRxSeconds(tRx_i, find(~isnan(gnssMeas.tRxSeconds(tRx_i, :)), 1)) <= 0
    for i = 1:GT_row
        if abs(GroundTruth(i,3)-gnssMeas.tRxSeconds(tRx_i, find(~isnan(gnssMeas.tRxSeconds(tRx_i, :)), 1))) <= 0.5
            % Align the first epoch
            break;
        end
    end
    % The ground truth file is 10 Hz
    GroundTruth_cache = GroundTruth(i:10:GT_row,:);  
%     % The ground truth file is 1 Hz
%     GroundTruth_cache = GroundTruth(i:GT_row,:);  
end

% Ground truth should not be longer than raw measurements  
[GT_row_cache, GT_clo_cache] = size(GroundTruth_cache);

if GT_row_cache > tRx_row - tRx_i +1 
    % Ground truth should not be longer than raw measurements 
    GT_row_cache = tRx_row - tRx_i +1;
elseif GT_row_cache < tRx_row - tRx_i +1
    % Ground truth should not be shorter than raw measurements 
    tRx_row = GT_row_cache + tRx_i - 1;
end
GroundTruth0 = GroundTruth_cache(1:GT_row_cache,:);


% GT_index = zeros(GT_row_cache,1);
% j = 0;
% % Find the corresponding ground truth
% % The loop of measurements
% for i = tRx_i:tRx_row    
%     % The loop of ground truth
%     for j = 1: GT_row_cache
%         if abs(GroundTruth_cache(j,3)-gnssMeas.tRxSeconds(i, find(~isnan(gnssMeas.tRxSeconds(i, :)), 1))) <= 0.5
%             GT_index(j) = 1;
%             break
%         end
%     end
% end
% GroundTruth0 = GroundTruth_cache(logical(GT_index),:);

[GT_row0, GT_clo0] = size(GroundTruth0);


%% Compute the pseudorange errors and the geometry matrix for smartphones
[PrmErr, Wls_dBc, h_end] = ComputePrErr(gpsPvt, GroundTruth0(:,7:9), tRx_i);

%% Estimate WLS estimation error of user clock bias
Hat_Wls_dBc = EstimateWlsDbc(gpsPvt, PrmErr);

%% Estimate the heading of the smartphone
Hat_Heading = EstimateHeading(gpsPvt);

%% Create Training Dataset file
f_3DSvPVTGT_test = fopen([dirName,'/SvPVT3D_Error_label_dynamic_',filename_date,'.csv'],'w');
fprintf(f_3DSvPVTGT_test,'%s , %s , %s , %s , %s , %s , %s , %s , %s, %s , %s , %s , %s, %s , %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n',...
'Epoch', 'PRN', 'SvX', 'SvY', 'SvZ','Dtsv','Elevation','Atmospheric Correction','CN0','R','WlsX','WlsY','WlsZ','Wlsdtu','GTX', 'GTY', 'GTZ',...
'Pseudorange Error Uncertainty','Pseudorange Error','DeltaDtu','Hat_DeltaDtu','Pseudorange Rate','WlsLon Degree','WlsLon Minute','WlsLon Second',...
'WlsLat Degree','WlsLat Minute','WlsLat Second', 'Unit Geometry Matrix N', 'Unit Geometry Matrix E', 'Unit Geometry Vector D', 'Pseudorange Error plus DeltaDtu','Azimuth', 'Pseudorange Residuals',...
'Smoothed Pseudorange Residuals','Heading of Smartphone N', 'Heading of Smartphone E', 'Heading of Smartphone D','Item of h','GtLon Degree','GtLon Minute','GtLon Second',...
'GtLat Degree','GtLat Minute','GtLat Second','EkfLon Degree','EkfLon Minute','EkfLon Second','EkfLat Degree','EkfLat Minute','EkfLat Second','normGtLon','normGtLat','normEkfLon','normEkfLat');     
NumSvPosR = size(gpsPvt.SvPosR, 1);
NumSvPosRS = size(gpsPvt.SvPosRS, 1);
if NumSvPosR ~= NumSvPosRS
    error("The number of EKF results don't equal to that of RTS results!");
end

% The size of the area
map_ll_0 =[-121.902, 37.325];
map_ll_m = [-121.88, 37.34];

for m = 1:NumSvPosR   
    n = gpsPvt.SvPosR(m, 1);
    if n < tRx_i
        % There is no ground truth for the measurements before tRx_i
        continue;
    end
    
    % ignore the first 120 data
    if n < 121 %|| n > gpsPvt.SvPosR(end, 1) - 100
%     if n < 1160 || n > 1275
        continue;
    end
        
    if isnan(gpsPvt.allXyzMMM(n,1)) || isnan(gpsPvt.allXyzMMM(n,2)) || isnan(gpsPvt.allXyzMMM(n,3))
        WlsX = 0;
        WlsY = 0;
        WlsZ = 0;
        Wlsdtu = 0;
        WlsLon = 0;
        WlsLat = 0;
    else
        WlsX = gpsPvt.allXyzMMM(n,1);
        WlsY = gpsPvt.allXyzMMM(n,2);
        WlsZ = gpsPvt.allXyzMMM(n,3);
        Wlsdtu = gpsPvt.allBcMeters(n);
        WlsLon = gpsPvt.allLlaDegDegM(n,2);
        WlsLat = gpsPvt.allLlaDegDegM(n,1);
    end
    
    % WLS Reuslts
    WlsLon_Degree = fix(WlsLon);
    WlsLon_Minute = fix((WlsLon-WlsLon_Degree)*60);
    WlsLon_Second = ((WlsLon-WlsLon_Degree)*60-WlsLon_Minute)*60;
    
    WlsLat_Degree = fix(WlsLat);
    WlsLat_Minute = fix((WlsLat-WlsLat_Degree)*60);
    WlsLat_Second = ((WlsLat-WlsLat_Degree)*60-WlsLat_Minute)*60;

    % Normalize GT positions: Longitude, Latitude, Altitude
    GroundTruthLla = [GroundTruth0(n-tRx_i+1,4), GroundTruth0(n-tRx_i+1,5), GroundTruth0(n-tRx_i+1,6)];
    
    normGroundTruthLl = (GroundTruthLla(1:2)-map_ll_0)./(map_ll_m-map_ll_0);
% 
%     % Normalize EKF positions: Longitude, Latitude, Altitude
%     if isnan(gpsPvt.allXyzMMMEKF(n,1)) || isnan(gpsPvt.allXyzMMMEKF(n,2)) || isnan(gpsPvt.allXyzMMMEKF(n,3))
%         EkfLon = 0;
%         EkfLat = 0;
%     else
%         EkfLon = gpsPvt.allLlaDegDegMEKF(n,2);
%         EkfLat = gpsPvt.allLlaDegDegMEKF(n,1);
%     end
%     normEkfLl = ([EkfLon, EkfLat]-map_ll_0)./(map_ll_m-map_ll_0);

    
    % Ground truth longitude and latitude 
    GtLon = GroundTruth0(n-tRx_i+1,4);
    GtLat = GroundTruth0(n-tRx_i+1,5);

    GtLon_Degree = fix(GtLon);
    GtLon_Minute = fix((GtLon-GtLon_Degree)*60);
    GtLon_Second = ((GtLon-GtLon_Degree)*60-GtLon_Minute)*60;
    
    GtLat_Degree = fix(GtLat);
    GtLat_Minute = fix((GtLat-GtLat_Degree)*60);
    GtLat_Second = ((GtLat-GtLat_Degree)*60-GtLat_Minute)*60;

    % EKF longitude and latitude
    if isnan(gpsPvt.allXyzMMMEKF(n,1)) || isnan(gpsPvt.allXyzMMMEKF(n,2)) || isnan(gpsPvt.allXyzMMMEKF(n,3))
        EkfLon = 0;
        EkfLat = 0;
    else
        EkfLon = gpsPvt.allLlaDegDegMEKF(n,2);
        EkfLat = gpsPvt.allLlaDegDegMEKF(n,1);
    end

    EkfLon_Degree = fix(EkfLon);
    EkfLon_Minute = fix((EkfLon-EkfLon_Degree)*60);
    EkfLon_Second = ((EkfLon-EkfLon_Degree)*60-EkfLon_Minute)*60;
    
    EkfLat_Degree = fix(EkfLat);
    EkfLat_Minute = fix((EkfLat-EkfLat_Degree)*60);
    EkfLat_Second = ((EkfLat-EkfLat_Degree)*60-EkfLat_Minute)*60;

    normEkfLl = ([EkfLon, EkfLat]-map_ll_0)./(map_ll_m-map_ll_0);
    
    % Groundtruth
    GroundTruthXyz = [GroundTruth0(n-tRx_i+1,7), GroundTruth0(n-tRx_i+1,8), GroundTruth0(n-tRx_i+1,9)];
    
    % Unit Geometry Vector using Wls results
    % 1. Calculate line of sight vectors from satellite to xo
    v = [WlsX, WlsY,WlsZ] - gpsPvt.SvPosR(m,4:6);%v(:,i) = vector from sv(i) to xyz0
    % 2. Calculate ranges from satellite to xo
    range = norm(v);
    % 3. Normalization
    if range == 0
        UnitGVector = [0, 0, 0];
    else
        UnitGVector = v./range; % line of sight unit vectors from sv to xo
    end
    % 4. Convert UnitGVector from ECEF to NED
    RE2N = RotEcef2Ned(WlsLat, WlsLon);
    UnitGVector = (RE2N*UnitGVector')';
        
    %Calculate pseudorange residuals for visible satellites
    % gpsPvt.SvPosR(m,7) * GpsConstants.LIGHTSPEED is the satellite clock bias
    
    PsResidual = PrmErr(m);
    h_end_item = h_end(m);
    PrMError_plus_bcError = gpsPvt.SvPosR(m,11) - norm(gpsPvt.SvPosR(m,4:6)-GroundTruthXyz)-Wlsdtu;
    PrMResi = gpsPvt.SvPosR(m,17) - norm(gpsPvt.SvPosR(m,4:6)-GroundTruthXyz);
    PrMResiS = gpsPvt.SvPosRS(NumSvPosRS + 1 - m) - norm(gpsPvt.SvPosR(m,4:6)-GroundTruthXyz);
%     PsResidual = norm(gpsPvt.SvPosR(m,4:6)-GroundTruthXyz)-gpsPvt.SvPosR(m,7) * GpsConstants.LIGHTSPEED;%gpsPvt.SvPosR(m,end-2)-norm(gpsPvt.SvPosR(m,4:6)-GroundTruthXyz)-Wlsdtu; % epsilon + delta_bc
    
    
    A = [gpsPvt.SvPosR(m,1), gpsPvt.SvPosR(m,3:6), gpsPvt.SvPosR(m,7)* GpsConstants.LIGHTSPEED, gpsPvt.SvPosR(m,8),gpsPvt.SvPosR(m,12)+gpsPvt.SvPosR(m,13),...
        gpsPvt.SvPosR(m,10),gpsPvt.SvPosR(m,14), WlsX, WlsY, WlsZ, Wlsdtu,GroundTruthXyz,gpsPvt.SvPosR(m,15),PsResidual,Wls_dBc(n), Hat_Wls_dBc(n), gpsPvt.SvPosR(m,16),...
        WlsLon_Degree, WlsLon_Minute, WlsLon_Second, WlsLat_Degree, WlsLat_Minute, WlsLat_Second, UnitGVector,PrMError_plus_bcError,gpsPvt.SvPosR(m,9), PrMResi,PrMResiS, Hat_Heading(n,:), h_end_item,...
        GtLon_Degree, GtLon_Minute, GtLon_Second, GtLat_Degree, GtLat_Minute, GtLat_Second,...
        EkfLon_Degree, EkfLon_Minute, EkfLon_Second, EkfLat_Degree, EkfLat_Minute, EkfLat_Second,...
        normGroundTruthLl, normEkfLl];
    
    fprintf(f_3DSvPVTGT_test,'%d , %d , %f , %f , %f , %f, %f , %f , %f , %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n',A');    
end
% fclose(f_3DSvPVTGT_test);

% Create data file for analysis
f_3DSvPVTGT_test = fopen([dirName,'/SvPVT3D_Error_label_dynamic_',filename_date,'.csv'],'r');
f_3DSvPVTGT_test_data = fopen([dirName,'/SvPVT3D_Error_label_dynamic_data_',filename_date,'.csv'],'w');
% Read the header
line = fgetl(f_3DSvPVTGT_test);
while ischar(line)
    
    if contains(line,'Epoch')   
        line = fgetl(f_3DSvPVTGT_test);
    else
        fprintf(f_3DSvPVTGT_test_data,'%s\n',line);
        line = fgetl(f_3DSvPVTGT_test);
    end
    
end
fclose(f_3DSvPVTGT_test);
fclose(f_3DSvPVTGT_test_data);



% %% Plot ADR
% [n,m] = size(gnssMeas.PrSM);
% 
% for i = 1:m
%    
%     figure;
%     for j = 2:n       
%         dADR(j-1) = gnssMeas.AdrM(j,i) - gnssMeas.AdrM(j-1,i);    
%     end
%     plot(dADR);
%     legend("PRN "+string(gnssMeas.Svid));
%     
% end

%% Compare WLS, EKF, MHE, Google's GNSS Analysis Software and Groundtruth
figure1 = figure;
% WLS
plot(gpsPvt.allLlaDegDegM(tRx_i:tRx_row-2,2),gpsPvt.allLlaDegDegM(tRx_i:tRx_row-2,1),'m','linewidth',4);hold on;
% MHE
plot(gpsPvt.allLlaDegDegMMHE(tRx_i:tRx_row-2,2),gpsPvt.allLlaDegDegMMHE(tRx_i:tRx_row-2,1),'r','linewidth',4);hold on;
% EKF
plot(gpsPvt.allLlaDegDegMEKF(tRx_i:tRx_row-2,2),gpsPvt.allLlaDegDegMEKF(tRx_i:tRx_row-2,1),'b','linewidth',4);hold on;
% RTS
plot(gpsPvt.allLlaDegDegMEKFS(tRx_i:tRx_row-2,2),gpsPvt.allLlaDegDegMEKFS(tRx_i:tRx_row-2,1),'c','linewidth',4);hold on;
% Google's GNSS Analysis Software
% plot(gnssDerived.WLS_lg(:),gnssDerived.WLS_la(:),'c','linewidth',2);hold on;
% plot(gnssDerived.ADR_lg(:),gnssDerived.ADR_la(:),'c','linewidth',2);hold on;
% Groundtruth
plot(GroundTruth0(1:GT_row0,4),GroundTruth0(1:GT_row0,5),'color','#008744','linewidth',4);hold on;


% figure1 = figure;
% % WLS
% plot(gpsPvt.allLlaDegDegM(1540:1590,2),gpsPvt.allLlaDegDegM(1540:1590,1),'m','linewidth',4);hold on;
% % MHE
% plot(gpsPvt.allLlaDegDegMMHE(1540:1590,2),gpsPvt.allLlaDegDegMMHE(1540:1590,1),'r','linewidth',4);hold on;
% % EKF
% plot(gpsPvt.allLlaDegDegMEKF(1540:1590,2),gpsPvt.allLlaDegDegMEKF(1540:1590,1),'b','linewidth',4);hold on;
% % RTS
% plot(gpsPvt.allLlaDegDegMEKFS(1540:1590,2),gpsPvt.allLlaDegDegMEKFS(1540:1590,1),'c','linewidth',4);hold on;
% % Google's GNSS Analysis Software
% % plot(gnssDerived.WLS_lg(:),gnssDerived.WLS_la(:),'c','linewidth',2);hold on;
% % plot(gnssDerived.ADR_lg(:),gnssDerived.ADR_la(:),'c','linewidth',2);hold on;
% % Groundtruth
% plot(GroundTruth0(1540:1590,4),GroundTruth0(1540:1590,5),'color','#008744','linewidth',4);hold on;


% figure1 = figure;
% % WLS
% plot(gpsPvt.allLlaDegDegM(1665:1690,2),gpsPvt.allLlaDegDegM(1665:1690,1),'m','linewidth',4);hold on;
% % MHE
% plot(gpsPvt.allLlaDegDegMMHE(1665:1690,2),gpsPvt.allLlaDegDegMMHE(1665:1690,1),'r','linewidth',4);hold on;
% % EKF
% plot(gpsPvt.allLlaDegDegMEKF(1665:1690,2),gpsPvt.allLlaDegDegMEKF(1665:1690,1),'b','linewidth',4);hold on;
% % RTS
% plot(gpsPvt.allLlaDegDegMEKFS(1665:1690,2),gpsPvt.allLlaDegDegMEKFS(1665:1690,1),'c','linewidth',4);hold on;
% % Google's GNSS Analysis Software
% % plot(gnssDerived.WLS_lg(:),gnssDerived.WLS_la(:),'c','linewidth',2);hold on;
% % plot(gnssDerived.ADR_lg(:),gnssDerived.ADR_la(:),'c','linewidth',2);hold on;
% % Groundtruth
% plot(GroundTruth0(1665:1690,4),GroundTruth0(1665:1690,5),'color','#008744','linewidth',4);hold on;

% figure1 = figure;
% % WLS
% plot(gpsPvt.allLlaDegDegM(1566:1583,2),gpsPvt.allLlaDegDegM(1566:1583,1),'m','linewidth',4);hold on;
% % MHE
% plot(gpsPvt.allLlaDegDegMMHE(1566:1583,2),gpsPvt.allLlaDegDegMMHE(1566:1583,1),'r','linewidth',4);hold on;
% % EKF
% plot(gpsPvt.allLlaDegDegMEKF(1566:1583,2),gpsPvt.allLlaDegDegMEKF(1566:1583,1),'b','linewidth',4);hold on;
% % RTS
% plot(gpsPvt.allLlaDegDegMEKFS(1566:1583,2),gpsPvt.allLlaDegDegMEKFS(1566:1583,1),'c','linewidth',4);hold on;
% % Google's GNSS Analysis Software
% % plot(gnssDerived.WLS_lg(:),gnssDerived.WLS_la(:),'c','linewidth',2);hold on;
% % plot(gnssDerived.ADR_lg(:),gnssDerived.ADR_la(:),'c','linewidth',2);hold on;
% % Groundtruth
% plot(GroundTruth0(1566:1583,4),GroundTruth0(1566:1583,5),'color','#008744','linewidth',4);hold on;



% legend("WLS","EKF",'MHE',"Google's GNSS Analysis Software Result","U-blox Receiver Ground Truth");
legend("WLS",'MHE',"EKF","RTS smoother","U-blox Receiver Ground Truth",'FontSize',20);
% legend("WLS","EKF","Google's GNSS Analysis Software Result","U-blox Receiver Ground Truth");
xlabel("Longitude",'FontSize',20);
ylabel("Latitude",'FontSize',20);
axis tight;
hold on;

% % SFO area
% I = imread([dirName '/SFO.png']);
% h = image('XData',[-122.35, -122.1],'YData',[37.5,37.4],'CData',I);%note the latitude (y-axis) is flipped in vertical direction % SFO


% SJC area
I = imread([dirName '/SJC.png']); 
axis([-121.902, -121.88, 37.325, 37.34]);
h = image('XData',[-121.902, -121.88],'YData',[37.34,37.325],'CData',I);%note the latitude (y-axis) is flipped in vertical direction

% % MTV area
% I = imread([dirName '/MTV.png']); 
% axis([-122.5, -121.8, 37.3, 37.7]);
% h = image('XData',[-122.5, -121.8],'YData',[37.7,37.3],'CData',I);%note the latitude (y-axis) is flipped in vertical direction
% % h = image('XData',[-122.5, -122.0],'YData',[37.7,37.4],'CData',I);%note the latitude (y-axis) is flipped in vertical direction

% % MTV-Generation
% I = imread([dirName '/MTVG.png']);hold on;
% axis([-122.10708,-122.07489,37.38197,37.40046]);  
% h = image('XData',[-122.10708,-122.07489],'YData',[37.40046,37.38197],'CData',I);%note the latitude (y-axis) is flipped in vertical direction

uistack(h,'bottom'); %move the image to the bottom of current stack
saveas(figure1,[dirName,'/Results.fig']);


%% Calculate the positioning errors

for n = tRx_i:tRx_row
   % WLS
    GroundTruthXyz = Lla2Xyz([GroundTruth0(n-tRx_i+1,5), GroundTruth0(n-tRx_i+1,4), GroundTruth0(n-tRx_i+1,6)]);
    delta_x(n,1) = gpsPvt.allXyzMMM(n,1)-GroundTruthXyz(1); % Delta X
    delta_y(n,1) = gpsPvt.allXyzMMM(n,2)-GroundTruthXyz(2); % Delta Y
    delta_z(n,1) = gpsPvt.allXyzMMM(n,3)-GroundTruthXyz(3); % Delta Z
    
    delta_xEKF(n,1) = gpsPvt.allXyzMMMEKF(n,1)-GroundTruthXyz(1); % Delta X
    delta_yEKF(n,1) = gpsPvt.allXyzMMMEKF(n,2)-GroundTruthXyz(2); % Delta Y
    delta_zEKF(n,1) = gpsPvt.allXyzMMMEKF(n,3)-GroundTruthXyz(3); % Delta Z
    
    delta_xEKFS(n,1) = gpsPvt.allXyzMMMEKFS(n,1)-GroundTruthXyz(1); % Delta X
    delta_yEKFS(n,1) = gpsPvt.allXyzMMMEKFS(n,2)-GroundTruthXyz(2); % Delta Y
    delta_zEKFS(n,1) = gpsPvt.allXyzMMMEKFS(n,3)-GroundTruthXyz(3); % Delta Z
    
    delta_xMHE(n,1) = gpsPvt.allXyzMMMMHE(n,1)-GroundTruthXyz(1); % Delta X
    delta_yMHE(n,1) = gpsPvt.allXyzMMMMHE(n,2)-GroundTruthXyz(2); % Delta Y
    delta_zMHE(n,1) = gpsPvt.allXyzMMMMHE(n,3)-GroundTruthXyz(3); % Delta Z
    
    % NED System
    %Calculate the rotation matrix to convert an ECEF vector to
    % North, East, Down coordinates, and vice-versa   
    RE2N = RotEcef2Ned(GroundTruth0(n-tRx_i+1,5), GroundTruth0(n-tRx_i+1,4));
    delta_v_ned = RE2N*[delta_x(n,1);delta_y(n,1);delta_z(n,1)];
    delta_n(n,1) = delta_v_ned(1);
    delta_e(n,1) = delta_v_ned(2);
    delta_d(n,1) = delta_v_ned(3);
    
    delta_v_ned_EKF = RE2N*[delta_xEKF(n,1);delta_yEKF(n,1);delta_zEKF(n,1)];
    delta_nEKF(n,1) = delta_v_ned_EKF(1);
    delta_eEKF(n,1) = delta_v_ned_EKF(2);
    delta_dEKF(n,1) = delta_v_ned_EKF(3);
    
    delta_v_ned_EKFS = RE2N*[delta_xEKFS(n,1);delta_yEKFS(n,1);delta_zEKFS(n,1)];
    delta_nEKFS(n,1) = delta_v_ned_EKFS(1);
    delta_eEKFS(n,1) = delta_v_ned_EKFS(2);
    delta_dEKFS(n,1) = delta_v_ned_EKFS(3);
    
    delta_v_ned_MHE = RE2N*[delta_xMHE(n,1);delta_yMHE(n,1);delta_zMHE(n,1)];
    delta_nMHE(n,1) = delta_v_ned_MHE(1);
    delta_eMHE(n,1) = delta_v_ned_MHE(2);
    delta_dMHE(n,1) = delta_v_ned_MHE(3);
    
    % Geodetic System
    GroundTruth = [GroundTruth0(n-tRx_i+1,5), GroundTruth0(n-tRx_i+1,4), GroundTruth0(n-tRx_i+1,6)];
    delta_lat(n,1) = gpsPvt.allLlaDegDegM(n,1)-GroundTruth(1); % Delta latitude
    delta_lon(n,1) = gpsPvt.allLlaDegDegM(n,2)-GroundTruth(2); % Delta longitude
    delta_alt(n,1) = gpsPvt.allLlaDegDegM(n,3)-GroundTruth(3); % Delta altitude
    
    delta_lat_EKF(n,1) = gpsPvt.allLlaDegDegMEKF(n,1)-GroundTruth(1); % Delta latitude
    delta_lon_EKF(n,1) = gpsPvt.allLlaDegDegMEKF(n,2)-GroundTruth(2); % Delta longitude
    delta_alt_EKF(n,1) = gpsPvt.allLlaDegDegMEKF(n,3)-GroundTruth(3); % Delta altitude
    
    delta_lat_EKFS(n,1) = gpsPvt.allLlaDegDegMEKFS(n,1)-GroundTruth(1); % Delta latitude
    delta_lon_EKFS(n,1) = gpsPvt.allLlaDegDegMEKFS(n,2)-GroundTruth(2); % Delta longitude
    delta_alt_EKFS(n,1) = gpsPvt.allLlaDegDegMEKFS(n,3)-GroundTruth(3); % Delta altitude
    
    delta_lat_MHE(n,1) = gpsPvt.allLlaDegDegMMHE(n,1)-GroundTruth(1); % Delta latitude
    delta_lon_MHE(n,1) = gpsPvt.allLlaDegDegMMHE(n,2)-GroundTruth(2); % Delta longitude
    delta_alt_MHE(n,1) = gpsPvt.allLlaDegDegMMHE(n,3)-GroundTruth(3); % Delta altitude
    
    
    
end

% % At the equator for longitude and for latitude: 1 degree = 111 km 
% delta_lat_m = delta_lat*111000;
% delta_long_m = delta_long*111000;
figure;
subplot(3,1,1);
plot(delta_x,'m','linewidth',2);hold on;
plot(delta_xEKF,'b','linewidth',2);hold on;
plot(delta_xMHE,'r','linewidth',2); hold on;
plot(delta_xEKFS,'c','linewidth',2);
axis([0, 1800, -100, 50]);
save('error_xWLS.mat','delta_x');
save('error_xEKF.mat','delta_xEKF');
save('error_xMHE.mat','delta_xMHE');
% save('error_xWLS_ND.mat','delta_x');
% save('error_xEKF_ND.mat','delta_xEKF');
% save('error_xMHE_ND.mat','delta_xMHE');
xlabel("Epoch (s)",'linewidth',2);
ylabel("Error on X axis (m)",'linewidth',2);
legend("WLS Results","EKF Results","MHE-WLS Results","RTS Smoothing Results");
% Mean_delta_lat_m = mean(delta_lat_m(40:130));
% Mean_delta_lat_m = mean(delta_lat_m,'omitnan');
% RMSE_lat_m = rms(delta_lat_m,'omitnan');
subplot(3,1,2);
plot(delta_y,'m','linewidth',2);hold on
plot(delta_yEKF,'b','linewidth',2);hold on;
plot(delta_yMHE,'r','linewidth',2);hold on;
plot(delta_yEKFS,'c','linewidth',2);
axis([0, 1800, -133, 50]);
save('error_yWLS.mat','delta_y');
save('error_yEKF.mat','delta_yEKF');
save('error_yMHE.mat','delta_yMHE');
% save('error_yWLS_ND.mat','delta_y');
% save('error_yEKF_ND.mat','delta_yEKF');
% save('error_yMHE_ND.mat','delta_yMHE');
xlabel("Epoch (s)",'linewidth',2);
ylabel("Error on Y axis  (m)",'linewidth',2);

subplot(3,1,3);
plot(delta_z,'m','linewidth',2);hold on
plot(delta_zEKF,'b','linewidth',2);hold on;
plot(delta_zMHE,'r','linewidth',2);hold on;
plot(delta_zEKFS,'c','linewidth',2);
axis([0, 1800, -40, 148]);
save('error_zWLS.mat','delta_z');
save('error_zEKF.mat','delta_zEKF');
save('error_zMHE.mat','delta_zMHE');
% save('error_zWLS_ND.mat','delta_z');
% save('error_zEKF_ND.mat','delta_zEKF');
% save('error_zMHE_ND.mat','delta_zMHE');
xlabel("Epoch (s)",'linewidth',2);
ylabel("Error on Z axis  (m)",'linewidth',2);
% Mean_delta_long_m = mean(delta_long_m(40:130));
% Mean_delta_long_m = mean(delta_long_m,'omitnan');
% RMSE_long_m = rms(delta_long_m,'omitnan');
% figure;
% subplot(2,1,1);
% plot(delta_lat,'g','linewidth',2);
% xlabel("Epoch (s)",'linewidth',2);
% ylabel("Error of Latitude (degree)",'linewidth',2);
% 
% subplot(2,1,2);
% plot(delta_long,'r','linewidth',2);
% xlabel("Epoch (s)",'linewidth',2);
% ylabel("Error of Longitude (degree)",'linewidth',2);


% Plot comparison of NED results
figure_NED = figure;
subplot(3,1,1);
plot(delta_n,'m','linewidth',2);hold on;
plot(delta_nMHE,'r','linewidth',2);hold on;
plot(delta_nEKF,'b','linewidth',2);hold on;
plot(delta_nEKFS,'c','linewidth',2);
xlabel("Epoch (s)",'FontSize',14);
ylabel("Error of North (m)",'FontSize',14);
legend("WLS","MHE-WLS","EKF","RTS Smoother",'FontSize',14);
Mean_delta_n_d = mean(delta_n,'omitnan');
RMSE_n = rms(delta_n,'omitnan');
RMSE_n_MHE = rms(delta_nMHE,'omitnan');
RMSE_n_EKF = rms(delta_nEKF,'omitnan');
RMSE_n_EKFS = rms(delta_nEKFS,'omitnan');


subplot(3,1,2);
plot(delta_e,'m','linewidth',2);hold on;
plot(delta_eEKF,'b','linewidth',2);hold on;
plot(delta_eMHE,'r','linewidth',2);hold on;
plot(delta_eEKFS,'c','linewidth',2);
xlabel("Epoch (s)",'FontSize',14);
ylabel("Error of East (m)",'FontSize',14);
% legend("WLS Results","EKF Results","MHE-WLS Results");
Mean_delta_e_d = mean(delta_e,'omitnan');
% Mean_delta_lon_EKF_d = mean(delta_lon_EKF,'omitnan');
% Mean_delta_lon_MHE_d = mean(delta_lon_MHE,'omitnan');
RMSE_e = rms(delta_e,'omitnan');
RMSE_e_MHE = rms(delta_eMHE,'omitnan');
RMSE_e_EKF = rms(delta_eEKF,'omitnan');
RMSE_e_EKFS = rms(delta_eEKFS,'omitnan');

subplot(3,1,3);
plot(delta_d,'m','linewidth',2);hold on;
plot(delta_dEKF,'b','linewidth',2);hold on;
plot(delta_dMHE,'r','linewidth',2);hold on;
plot(delta_dEKFS,'c','linewidth',2);
xlabel("Epoch (s)",'FontSize',14);
ylabel("Error of Down (Meter)",'FontSize',14);
% legend("WLS Results","EKF Results","MHE-WLS Results");
Mean_delta_d_d = mean(delta_d,'omitnan');
% Mean_delta_alt_EKF_d = mean(delta_alt_EKF,'omitnan');
% Mean_delta_alt_MHE_d = mean(delta_alt_MHE,'omitnan');
RMSE_d = rms(delta_d,'omitnan');
RMSE_d_MHE = rms(delta_dMHE,'omitnan');
RMSE_d_EKF = rms(delta_dEKF,'omitnan');
RMSE_d_EKFS = rms(delta_dEKFS,'omitnan');
saveas(figure_NED,[dirName,'/NED_Results.fig']);

% Plot comparison of LLA results
figure;
subplot(3,1,1);
plot(delta_lat,'m','linewidth',2);hold on;
plot(delta_lat_EKF,'b','linewidth',2);hold on;
plot(delta_lat_MHE,'r','linewidth',2);hold on;
plot(delta_lat_EKFS,'c','linewidth',2);
xlabel("Epoch (s)",'linewidth',2);
ylabel("Error of Latitude (Degree)",'linewidth',2);
legend("WLS Results","EKF Results","MHE-WLS Results","RTS Smoothing Results");
Mean_delta_lat_d = mean(delta_lat,'omitnan');
Mean_delta_lat_EFK_d = mean(delta_lat_EKF,'omitnan');
Mean_delta_lat_MHE_d = mean(delta_lat_MHE,'omitnan');
RMSE_lat_d = rms(delta_lat,'omitnan');
RMSE_lat_EKF_d = rms(delta_lat_EKF,'omitnan');
RMSE_lat_MHE_d = rms(delta_lat_MHE,'omitnan');

subplot(3,1,2);
plot(delta_lon,'m','linewidth',2);hold on;
plot(delta_lon_EKF,'b','linewidth',2);hold on;
plot(delta_lon_MHE,'r','linewidth',2);hold on;
plot(delta_lon_EKFS,'c','linewidth',2);
xlabel("Epoch (s)",'linewidth',2);
ylabel("Error of Longitude (Degree)",'linewidth',2);
Mean_delta_lon_d = mean(delta_lon,'omitnan');
Mean_delta_lon_EKF_d = mean(delta_lon_EKF,'omitnan');
Mean_delta_lon_MHE_d = mean(delta_lon_MHE,'omitnan');
RMSE_lon_d = rms(delta_lon,'omitnan');
RMSE_lon_EKF_d = rms(delta_lon_EKF,'omitnan');
RMSE_lon_MHE_d = rms(delta_lon_MHE,'omitnan');

subplot(3,1,3);
plot(delta_alt,'m','linewidth',2);hold on;
plot(delta_alt_EKF,'b','linewidth',2);hold on;
plot(delta_alt_MHE,'r','linewidth',2);hold on;
plot(delta_alt_EKFS,'c','linewidth',2);
xlabel("Epoch (s)",'linewidth',2);
ylabel("Error of Altitude (Meter)",'linewidth',2);
Mean_delta_alt_d = mean(delta_alt,'omitnan');
Mean_delta_alt_EKF_d = mean(delta_alt_EKF,'omitnan');
Mean_delta_alt_MHE_d = mean(delta_alt_MHE,'omitnan');
RMSE_alt_d = rms(delta_alt,'omitnan');
RMSE_alt_EKF_d = rms(delta_alt_EKF,'omitnan');
RMSE_alt_MHE_d = rms(delta_alt_MHE,'omitnan');

%% Write WLS results into a file
% offset = 120;
offset = 514;
tRx_row = 650;
tRx_i = tRx_i + offset;
f_3D_WLS_PVT = fopen([dirName,'/WLS_PVT_3D_test.csv'],'w');
index = (tRx_i:tRx_row)';
M = [index, gpsPvt.allLlaDegDegM(tRx_i:tRx_row,2),gpsPvt.allLlaDegDegM(tRx_i:tRx_row,1),gpsPvt.allLlaDegDegM(tRx_i:tRx_row,3)];
fprintf(f_3D_WLS_PVT,'%f , %f, %f, %f\n',M');
fclose(f_3D_WLS_PVT);


%% Write EKF results into a file
f_3D_EKF_PVT = fopen([dirName,'/EKF_PVT_3D_test.csv'],'w');
index = (tRx_i:tRx_row)';
M = [index, gpsPvt.allLlaDegDegMEKF(tRx_i:tRx_row,2),gpsPvt.allLlaDegDegMEKF(tRx_i:tRx_row,1),gpsPvt.allLlaDegDegMEKF(tRx_i:tRx_row,3)];
fprintf(f_3D_EKF_PVT,'%f , %f, %f, %f\n',M');
fclose(f_3D_EKF_PVT);

%% Write RTS Smoothing results into a file
f_3D_RTS_PVT = fopen([dirName,'/RTS_PVT_3D_test.csv'],'w');
index = (tRx_i:tRx_row)';
M = [index, gpsPvt.allLlaDegDegMEKFS(tRx_i:tRx_row,2),gpsPvt.allLlaDegDegMEKFS(tRx_i:tRx_row,1),gpsPvt.allLlaDegDegMEKFS(tRx_i:tRx_row,3)];
fprintf(f_3D_RTS_PVT,'%f , %f, %f, %f\n',M');
fclose(f_3D_RTS_PVT);

%% Write MHE results into a file
f_3D_MHE_PVT = fopen([dirName,'/MHE_PVT_3D_test.csv'],'w');
index = (tRx_i:tRx_row)';
M = [index, gpsPvt.allLlaDegDegMMHE(tRx_i:tRx_row,2),gpsPvt.allLlaDegDegMMHE(tRx_i:tRx_row,1),gpsPvt.allLlaDegDegMMHE(tRx_i:tRx_row,3)];
fprintf(f_3D_MHE_PVT,'%f , %f, %f, %f\n',M');
fclose(f_3D_MHE_PVT);

%% Write ground truth results into a file
f_3D_GT_PVT = fopen([dirName,'/GT_PVT_3D.csv'],'w');
GT_row0 = 650;
index = (1+offset:GT_row0)';
M = [index, GroundTruth0(1+offset:GT_row0,4),GroundTruth0(1+offset:GT_row0,5),GroundTruth0(1+offset:GT_row0,6)];
fprintf(f_3D_GT_PVT,'%f , %f, %f, %f\n',M');
fclose(f_3D_GT_PVT);


% %% Write WLS results and U-blox results into a file
% f_3D_WLS_uBlox_PVT = fopen('WLS_uBlox_PVT_3D.csv','w');
% M = [gpsPvt.allLlaDegDegM(tRx_i:tRx_row,2),gpsPvt.allLlaDegDegM(tRx_i:tRx_row,1),GroundTruth0(1:GT_row0,4),GroundTruth0(1:GT_row0,5)];
% fprintf(f_3D_WLS_uBlox_PVT,'%f , %f , %f , %f\n',M');
% fclose(f_3D_WLS_uBlox_PVT);

% %% plot Pvt results
% h4 = figure;
% ts = 'Raw Pseudoranges, Weighted Least Squares solution';
% PlotPvt(gpsPvt,prFileName,param.llaTrueDegDegM,ts); drawnow;
% h5 = figure;
% PlotPvtStates(gpsPvt,prFileName);
% 
% %% Plot Accumulated Delta Range 
% if any(any(isfinite(gnssMeas.AdrM) & gnssMeas.AdrM~=0))
%     [gnssMeas]= ProcessAdr(gnssMeas);
%     h6 = figure;
%     PlotAdr(gnssMeas,prFileName,colors);
%     [adrResid]= GpsAdrResiduals(gnssMeas,allGpsEph,param.llaTrueDegDegM);drawnow
%     h7 = figure;
%     PlotAdrResids(adrResid,gnssMeas,prFileName,colors);
% end
%% end of ProcessGnssMeasScript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2016 Google Inc.
% Copyright 2022 Nanyang Technological University
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
