function corrProj_7_6_overlapHistograms_SEAN_MOD(cellRegion)
% function corrProj_7_6_overlapHistograms(cellRegion)
%
% PURPOSE:
%  To plot the histograms of cross correlation coefficient
%  distributions by spatial overlap.
%
% INPUT:
%  cellRegion = uber data struct for each cell pair, outputted by corrProj_6...
%
% OUTPUT:
%  Figures
%
% JBT 8/2017
% Colgin Lab


stateNames = {'RUN', 'REM', 'NREM'};
groupNames = {'High', 'Mid', 'Far'};
% stateCols = {'Green', 'Gold', 'Purple'};
stateCols = {[0.0 0.5 0.0], [1.0 0.54688 0.0],[0.29297 0.0 0.50781]};
ratioCuts = [ 0.7 1.0;... %NEAR
    0.3 0.7;... %MID
    0.0 0.3]; %FAR

reg = 1; %only an MEC plot

for cp = 1:length(cellRegion(reg).cellPair)
    allSpatPhiMags(cp) = cellRegion(reg).cellPair(cp).state(1).relSpatPhiMag;
end
ratioCuts = ratioCuts .* max(allSpatPhiMags);
max(allSpatPhiMags)
gCntr = [1 1 1]; 
for cp = 1:length(cellRegion(reg).cellPair)
    cpSpatPhiMag = allSpatPhiMags(cp); 
%     orientRatio = cellRegion(reg).cellPair(cp).state(1).orientRatio;
    grpNum = [];
    if cpSpatPhiMag >= ratioCuts(1,1)% && orientRatio >= ratioCuts(1,1)
        grpNum = 1; %NEAR
    elseif cpSpatPhiMag >= ratioCuts(2,1) && cpSpatPhiMag <= ratioCuts(2,2)  %&& ...
%             orientRatio >= ratioCuts(2,1) && orientRatio <= ratioCuts(2,2)
        grpNum = 2; %MID
    elseif cpSpatPhiMag <= ratioCuts(3,1)% && orientRatio <= ratioCuts(3,1)
        grpNum = 3; %FAR
    end
    
    if ~isempty(grpNum)
        
        for s = 1:3
            midSum = cellRegion(reg).cellPair(cp).state(s).midSum(1);
            group(grpNum).xCorrs(gCntr(grpNum),s) = midSum; 
        end
        gCntr(grpNum) = gCntr(grpNum) + 1; 
    end
end

figure('Position', [ 377    65   608   606])
for s = 1:3
    for g = 1:3
        subplot(3,3,(s-1)*3+g);
        [cnts,edges] = histcounts(group(g).xCorrs(:,s), 'Normalization', 'Probability', 'BinWidth', .01, 'BinLimits', [0 .1]); 
        edges = edges(2:end)-mean(diff(edges)); 
        bg = bar(edges,cnts); 
%         set(bg, 'FaceColor', rgb(stateCols{s})); 
        set(bg, 'FaceColor', (stateCols{s})); 
        xlim([0 .1]); 
        ylim([0 1]); 
        if s == 1
            title([groupNames{g} ' Overlap']); 
        end
        if g == 1
            ylabel({stateNames{s} ; 'Probability'}); 
        end
        if s == 3 && g == 2
            xlabel('Cross-Correlation Coefficient'); 
        end
    end
end



end%fnctn