function PrmMWA = MWAver(gnssMeas,M,N,WindowSize)

PrmMWA = zeros(N,M)+NaN;

for indexSV = 1:M
    i = 1;
    while i<N+1
        %find the first non-NaN row
        while isnan(gnssMeas.PrM(i,indexSV))
            i = i+1;
            if i>N
                break;
            end
        end
        
        if i>N
            break;
        end
        
        PrmMWA(i,indexSV)  = gnssMeas.PrM(i,indexSV); % initialize the pseudorange
        
        j = i -1;
        if j <= 0
            i = i+1;
            if i>N
                break;
            else
                continue;
            end
        end
        % Determine whether the current pseudorange can be averaged
        while ~isnan(gnssMeas.PrM(j ,indexSV)) && (i-j) < WindowSize
            j = j - 1;
            if j <= 0
                break;
            end
        end
        
        if j <= 0
            i = i + 1;
            continue;
        end
        
        if (i-j) == WindowSize && ~isnan(gnssMeas.PrM(j ,indexSV))
            
            sum_dPrsM = 0;
            sum_PrsM = 0;
            
            % Average Delta-Pseudoranges
            for k = j : i-1
                sum_dPrsM = sum_dPrsM + (gnssMeas.PrM(k+1 ,indexSV) - gnssMeas.PrM(k ,indexSV));
            end
            
            for k = j : i
                sum_PrsM = sum_PrsM + gnssMeas.PrM(k ,indexSV) + (i-k)*sum_dPrsM / WindowSize;
            end
            
            % Average Pseudoranges
            
            PrmMWA(i,indexSV) = sum_PrsM / (WindowSize+1);
        end    
        i = i+1;           
    end
end

end