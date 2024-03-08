function [xHat,Sum_numSvs,xHat_k_N] = MHEstimatorTest(CurrentEpoch,gnssMeas,allGpsEph,WindowSize,xo,weekNum,GT_data)

xHat = []; % state updates at the kth epoch
xHat_k_N = []; % state updates at the k-Nth epoch

% MHE based WLS Positioning for the i+L epoch

% Determine the current position of the MHE window
if CurrentEpoch <= WindowSize
    L = CurrentEpoch-1;
else
    L = WindowSize;
end


% Allocate space to the measurements in a window
prs_stack = cell([L+1,1]);

% Allocate space to the Ephemeris in a window
gpsEph_stack = cell([L+1,1]);

iValid0 = (isfinite(gnssMeas.PrM(CurrentEpoch-L:CurrentEpoch,:)));

Sum_numSvs = 0;

for j=1:L+1
    
    iValid = find(iValid0(j,:));
    
    svid = gnssMeas.Svid(iValid)';
    
    [gpsEph,iSv] = ClosestGpsEph(allGpsEph,svid,gnssMeas.FctSeconds(CurrentEpoch-L+j-1));
    
    gpsEph_stack(j) = {gpsEph};
    
    svid = svid(iSv); %svid for which we have ephemeris
    
    numSvs = length(svid); %number of satellites this epoch
    
    Sum_numSvs = Sum_numSvs+numSvs;% Total number of visible sv in the window
    
%     gpsPvt.numSvs(CurrentEpoch-L+j-1) = numSvs;
    
    %for those svIds with valid ephemeris, pack prs matrix for WlsNav
    prM = gnssMeas.PrM(CurrentEpoch-L+j-1,iValid(iSv))';
    if length(prM) >= 4
        index_GT = GT_data(:,2) == (CurrentEpoch-L+j-1);
        prM = GT_data(index_GT,end);      
    end
    prSigmaM= gnssMeas.PrSigmaM(CurrentEpoch-L+j-1,iValid(iSv))';
    prrMps  = gnssMeas.PrrMps(CurrentEpoch-L+j-1,iValid(iSv))';
    prrSigmaMps = gnssMeas.PrrSigmaMps(CurrentEpoch-L+j-1,iValid(iSv))';
    tRx = [ones(numSvs,1)*weekNum(CurrentEpoch-L+j-1),gnssMeas.tRxSeconds(CurrentEpoch-L+j-1,iValid(iSv))'];
    
    prs = [tRx, svid, prM, prSigmaM, prrMps, prrSigmaMps];
    prs_stack(j) = {prs};
end

if Sum_numSvs < 4
    return;%skip to next epoch(next window)
end

% MHE-WLS
[xHat,~,~,H,~,xHat_k_N] = WlsPvtMHE2ndOrder(prs_stack,gpsEph_stack,xo);

end