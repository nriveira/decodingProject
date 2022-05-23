function fmr1CircTrack_x_thetaDecodeByPhase(group)
% function fmr1CircTrack_x_thetaDecodeByPhase(group)
%
% JB Trimper
% 01/2021
% Colgin Lab

% NOTE: This function will produce a summary plot with slightly less theta
% cycles than the DecodebyTime code. This is due to the way that
% get_theta_phase_times output is formatted.


plotEachCycle = 0; %set to 1 to plot each cycle (waves, spikes, and PPM)
%               ****   If you choose to plot each cycle, you'll need to
%         enter a manual stop because there are hundreds and your PC will crash  ****


figCntr = 0; %initialize
maxNumFigs = 10; %stop point

saveFigs = 0; %Set to 1 to save figures
saveDir =  'E:\FMR1_CIRCTRACK\RESULTS\updatedFigs_033121\decodeByPhase';

shiftPlots = 1; %set to 1 to have individual cycle plots shifted so rat's current location is always plotted at 180
%                        ***This needs to be set to 1 for the average plot to be interpretable ***
%              *** Only set it to zero for the sake of visual inspection of individual theta cycle plots ***


% Decoding parameters
sampRate = 12; % # of phases I've broken down each theta cycle into
bayesWin = 1/12; % 1 phase bin
bayesStep = 1/12;

minNumU = 3; %minimum # of cells that need to fire at least minNumSpks each
minNumSpks = 3; %these numbers from Zheng et al

minThetaAmp = -inf; %Minimum theta amplitude --> zscore((abs(hilbert(thetaFilteredLfp)))

runThresh = 5; %cm/s

radBinCtrs = group(2).rat(1).day(1).binCtrs; %doesn't change across days/rats
newRewLoc = 0; %used to get all plots to line up as though rewards are at 0 & 180
[~,newRewInd] = min(abs(circ_dist(deg2rad(radBinCtrs), deg2rad(newRewLoc))-0));

curDir = pwd; %Note where we are so we can return to it later

meanAmps = [];

for g = 1:2
    fprintf('Group %d\n', g);
    
    if g == 1
        r = 1;
    else
        r = 2;
    end
    %     for r = 2
    fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);
    
    
    
    dayNums = 1:length(group(g).rat(r).day);
    if g == 2
        dayNums = 2; %So we're only considering the KO rat's day w/ most cells
    end
    for d = dayNums
        %         for d = 1:length(group(g).rat(r).day);
        fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day));
        
        %For saving all PPMs, error, and slope across days
        dayPpm = [];
        decErr = [];
        slopes = [];
        
        % Get ratemaps
        try
            rateMaps = zeros(length(group(g).rat(r).day(d).xBeginUnitInfo), length(group(g).rat(r).day(d).xBeginUnitInfo(1).smRateMap));
        catch
            keyboard
        end
        
        badU = [];
        uIDs = zeros(length(group(g).rat(r).day(d).xBeginUnitInfo),2);
        for u = 1:length(group(g).rat(r).day(d).xBeginUnitInfo)
            if max(group(g).rat(r).day(d).xBeginUnitInfo(u).rateMap)>=1
                rateMaps(u,:) = group(g).rat(r).day(d).xBeginUnitInfo(u).smRateMap; %Smoothed ratemap
                uIDs(u,:) = group(g).rat(r).day(d).xBeginUnitInfo(u).ID;
            else
                badU = [badU u]; %#ok
            end
        end
        rateMaps(badU,:) = [];
        rateMaps(rateMaps==0) = 0.0001; %get rid of zeros because our Bayesian decoder can't handle 'em.
        uIDs(badU,:) = [];
        
        if size(rateMaps,1) > 1 %If there are at least 2 units...
            
            % Shift the rate maps so, across all data, reward locations are always plotted at 0 & 180
            %   So plots here will always look like rat goes from 0 to 360
            rewLoc = group(g).rat(r).day(d).rewLocs(1);
            [~,rewInd] = min(abs(circ_dist(deg2rad(radBinCtrs), deg2rad(rewLoc))-0));
            shiftVal = newRewInd - rewInd;
            %             rateMaps = circshift(rateMaps, [0 shiftVal]);
            
            
            
            allInds = [];
            goodPhiInds = [];
            goodRs = [];
            goodSpks = [];
            goodUnits = [];
            goodLfpInds = [];
            
            
            for b = 1:4
                fprintf('\t\t\tBegin %d\n', b)
                
                %Get rat's speed
                radPos = group(g).rat(r).day(d).begin(b).radPos;
                instRs = get_runspeed(group(g).rat(r).day(d).begin(b).coords);
                smRs = smooth_runspeed(instRs);
                
                cd(group(g).rat(r).day(d).begin(b).dir)
                
                tetNum = group(g).rat(r).day(d).thetaTet;
                %                 if g == 1 & r == 1
                %                    tetNum = 10;
                forFigName = 12;
                %                 end
                
                lfpRoot = ['CSC' num2str(tetNum)];
                lfpStruct = read_in_lfp(['CSC' num2str(tetNum) '.ncs']);
                
                load([lfpRoot '_broadThetaLfp.mat']); %#ok
                lfpStruct.broadThetaLfp = filtLfp;
                
                load([lfpRoot '_narrowThetaLfp.mat']);  %#ok
                lfpStruct.narrowThetaLfp = filtLfp;
                thetaFiltLfpAmp = zscore(abs(hilbert(filtLfp)));
                
                cd(curDir);
                
                
                % Upsample the phases so not just peak/fall/trough/rise
                phiTms = get_theta_phase_times(lfpStruct);
                
                morePhiTms = zeros(12,size(phiTms,2));
                morePhiTms(1,:) = phiTms(1,:);
                morePhiTms(4,:) = phiTms(2,:);
                morePhiTms(7,:) = phiTms(3,:);
                morePhiTms(10,:) = phiTms(4,:);
                
                phiDifs = (phiTms(2,:) - phiTms(1,:)) ./ 3;
                morePhiTms(2,:) = morePhiTms(1,:) + phiDifs;
                morePhiTms(3,:) = morePhiTms(1,:) + phiDifs.*2;
                
                phiDifs = (phiTms(3,:) - phiTms(2,:)) ./ 3;
                morePhiTms(5,:) = morePhiTms(4,:) + phiDifs;
                morePhiTms(6,:) = morePhiTms(4,:) + phiDifs.*2;
                
                phiDifs = (phiTms(4,:) - phiTms(3,:)) ./ 3;
                morePhiTms(8,:) = morePhiTms(7,:) + phiDifs;
                morePhiTms(9,:) = morePhiTms(7,:) + phiDifs.*2;
                
                phiDifs = (phiTms(1,2:end) - phiTms(4,1:end-1)) ./ 3;
                morePhiTms(11,1:end-1) = morePhiTms(10,1:end-1) + phiDifs;
                morePhiTms(11,end) = morePhiTms(10,end) + phiDifs(end);
                morePhiTms(12,1:end-1) = morePhiTms(10,1:end-1) + phiDifs.*2;
                morePhiTms(12,end) = morePhiTms(10,end) + phiDifs(end).*2;
                
                phiTms = morePhiTms;

                
                lapInds = group(g).rat(r).day(d).begin(b).lapInds;
                lapTms = group(g).rat(r).day(d).begin(b).lapTms;
                
                
                for lp = 1:size(lapInds,1)
                    
                    startLpTm = lapTms(lp,1); %Referenced to video
                    endLpTm = lapTms(lp,2);
                    
                    
                    % DEFINE EDGE OF CYCLE AS TROUGHS
                    lapPhiInds = find(phiTms(7,:)>=startLpTm & phiTms(7,:)<=endLpTm); %Use troughs
                    lapPhiTms = phiTms(:,lapPhiInds); %#ok
                    
                    allInds = [allInds 2:size(lapPhiTms,2)-1];
                    
                    for p = 2:size(lapPhiTms,2)-1
                        
                        startCycle = lapPhiTms(7,p);
                        endCycle = lapPhiTms(7,p+1);
                        
                        if p-1>0 & p+2<=size(lapPhiTms,2) %#ok
                            
                            goodPhiInds = [goodPhiInds p];
                            
                            startWin = lapPhiTms(7,p-1); %now covering 3 cycles
                            endWin = lapPhiTms(7,p+2);
                            winPhiTms = [lapPhiTms(7:end,p-1); lapPhiTms(:,p); lapPhiTms(:,p+1); lapPhiTms(1:6,p+2)];
                            
                            lfpInds = find(lfpStruct.ts>=startCycle & lfpStruct.ts<=endCycle);
                            meanAmp = mean(thetaFiltLfpAmp(lfpInds)); %#ok
                            meanAmps = [meanAmps meanAmp]; %#ok
                            
                            meanRs = mean(smRs(smRs(:,1)>=startCycle & smRs(:,1)<=endCycle,2));
                            if meanRs >= runThresh  &&  meanAmp >= minThetaAmp
                                
                                goodRs = [goodRs p];
                                
                                
                                %Pre-allocate for spike raster
                                spkRstr = zeros(size(uIDs,1), 36); %36 = 3 cycles with phase split into 12
                                
                                cycleSpkPhis = cell(1,size(uIDs,1));
                                winSpkPhis = cell(1,size(uIDs,1));
                                
                                uCntr = 0;
                                for u = 1:length(group(g).rat(r).day(d).begin(b).unit)
                                    uID = group(g).rat(r).day(d).begin(b).unit(u).ID;
                                    if ismember(uID, uIDs, 'row') %if the unit wasn't discarded due to low firing ratemap
                                        uCntr = uCntr + 1;
                                        
                                        allSpkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                                        tmpCycleSpkTms = allSpkTms(allSpkTms>=startCycle & allSpkTms<=endCycle);
                                        
                                        if length(tmpCycleSpkTms)>=minNumSpks %at least 3 spks
                                            
                                            goodSpks = [goodSpks p];
                                            
                                            cycleSpkTms = tmpCycleSpkTms;
                                            winSpkTms = allSpkTms(allSpkTms>=startWin & allSpkTms<=endWin);
                                            
                                            for st = 1:length(cycleSpkTms)
                                                tmpTm = cycleSpkTms(st);
                                                phiInd = find(winPhiTms<=tmpTm, 1, 'Last');
                                                cycleSpkPhis{uCntr}(st) = phiInd;
                                            end
                                            
                                            winSpkPhis{uCntr} = zeros(1,36); %initialize for this unit
                                            for st = 1:length(winSpkTms)
                                                tmpTm = winSpkTms(st);
                                                phiInd = find(winPhiTms<=tmpTm, 1, 'Last');
                                                
                                                %                                                 winSpkPhis{uCntr}(st) = phiInd;  Changed so mult spikes in same bin can be plotted below
                                                winSpkPhis{uCntr}(phiInd) = winSpkPhis{uCntr}(phiInd) + 1; %add this to num spikes that occured in this bin
                                                spkRstr(uCntr,phiInd) = spkRstr(uCntr,phiInd) + 1;
                                            end
                                            
                                        end
                                    end
                                    
                                end %unit
                                
                                
                                
                                if sum(sum(spkRstr,2)>=minNumSpks)>=minNumU %If # of cells in spkRstr with >=minNumSpks is >=minNumU, then decode
                                    
                                    goodUnits = [goodUnits p];
                                    
                                    %Plot the cycle
                                    lfpInds = find(lfpStruct.ts>=startWin & lfpStruct.ts<=endWin);
                                    if length(lfpInds)>750
                                        lfpInds(751:end) = [];
                                    end
                                    
                                    if lfpInds(1)>0   &&   lfpInds(end)<length(lfpStruct.ts)
                                        
                                        goodLfpInds = [goodLfpInds p];
                                        
                                        % Do the decoding
                                        ppm = BayesianDecoder(spkRstr,rateMaps,bayesWin, bayesStep, sampRate);
                                        ppm(isnan(ppm)) = 1/size(rateMaps,2);
                                        posInds = find(radPos(:,1)>=startWin & radPos(:,1)<=endWin);
                                        posInds = [posInds(1)-1; posInds; posInds(end)+1]; %#ok
                                        winMids = winPhiTms;
                                        
                                        % Get rat's position for every bayes window
                                        posMat = radPos(posInds,:);
                                        allTms = winMids;
                                        knownPts = posMat(:,1);
                                        knownPos = posMat(:,2);
                                        allPos = interp1(knownPts, knownPos, allTms);
                                        if isnan(allPos(end))
                                            allPos(end) = allPos(end-1);
                                        end
                                        
                                        % Shift PPM so always relative to rat's position at 180
                                        if shiftPlots == 1
                                            for c = 1:length(allTms)
                                                shiftDeg = 180-allPos(c);
                                                shiftVal = round(shiftDeg/mean(diff(group(g).rat(r).day(d).binCtrs))); % the divisor is the radial bin size in degrees
                                                ppm(:,c) = circshift(ppm(:,c), shiftVal);
                                            end
                                        end
                                        
                                        
                                        % Rat's most probable position given the PPM
                                        [maxVals,maxInds] = max(ppm);
                                        maxInds(maxVals==1/size(rateMaps,2)) = NaN;
                                        cycleInds = 13:24;
                                        indDist = maxInds(cycleInds) - 36.5; % 36.5 because that's where 180* would be and
                                        %                                           the rat's position is shift to always be at 180
                                        radDist = indDist(~isnan(indDist)) .* 5; % because that's the radial bin size
                                        decErr = [decErr std(radDist)]; %#ok
                                        radBinCtrs = group(g).rat(r).day(d).binCtrs; %doesn't change across days/rats
                                        knownX = find(~isnan(maxInds(cycleInds)));
                                        cyclePos = maxInds(cycleInds);
                                        cyclePos = radBinCtrs(cyclePos(~isnan(cyclePos)));
                                        coefs = polyfit(knownX, cyclePos, 1);
                                        fitLineY = polyval(coefs, 1:length(cycleInds));
                                        slopes = [slopes coefs(1)]; %#ok
                                        
                                        
                                        % Store across begins
                                        if isempty(dayPpm)
                                            dayPpm = ppm;
                                        else
                                            dayPpm = cat(3,dayPpm, ppm);
                                        end
                                        
                                        
                                        if plotEachCycle == 1 && figCntr < maxNumFigs
                                            figure('Position', [216   147   795   633]);
                                            figCntr = figCntr + 1;
                                            
                                            subplot(3,1,1);
                                            plot(lfpStruct.ts(lfpInds(1):lfpInds(end)), lfpStruct.data(lfpInds(1):lfpInds(end)))
                                            hold on;
                                            plot(lfpStruct.ts(lfpInds(1):lfpInds(end)), lfpStruct.broadThetaLfp(lfpInds(1):lfpInds(end)))
                                            plot(lfpStruct.ts(lfpInds(1):lfpInds(end)), lfpStruct.narrowThetaLfp(lfpInds(1):lfpInds(end)))
                                            yBnds = get(gca, 'YLim');
                                            line([startCycle startCycle], [yBnds(1) yBnds(2)], 'Color', [.5 .5 .5])
                                            line([endCycle endCycle], [yBnds(1) yBnds(2)], 'Color', [.5 .5 .5])
                                            xlim(lfpStruct.ts([lfpInds(1) lfpInds(end)]))
                                            ylabel('Amplitude (\mu V)');
                                            title([group(g).name ', Rat ' num2str(r) ', Day ' num2str(d) ', Begin' num2str(b) ', Lap ' num2str(lp) ', Cycle ' num2str(p)])
                                            fix_font;
                                            xlabel('Time (s)');
                                            
                                            subplot(3,1,2);
                                            hold on;
                                            uCntr = 0;
                                            for u = 1:length(cycleSpkPhis)
                                                if ~isempty(cycleSpkPhis{u})
                                                    uCntr = uCntr + 1;
                                                    %                                                     for st =
                                                    %                                                     1:length(winSpkPhis{u})
                                                    %                                                         keyboard
                                                    %                                                         line([winSpkPhis{u}(st) winSpkPhis{u}(st)], [uCntr-.3 uCntr+.3])
                                                    %                                                     end
                                                    spkBins = find(winSpkPhis{u}); %find bins that had spikes
                                                    for i = 1:length(spkBins)
                                                        curBin = spkBins(i); %get the spike bin we want to look at
                                                        numBinSpks = winSpkPhis{u}(curBin);
                                                        int = 1/numBinSpks;
                                                        xVals = curBin:int:curBin + 1;
                                                        for st = 1:numBinSpks
                                                            
                                                            line([xVals(st) xVals(st)], [uCntr-.3 uCntr+.3])
                                                        end %spike times in this bin
                                                    end %all phiInds that had spikes
                                                end
                                            end
                                            ylim([0.1 uCntr+.9]);
                                            yBnds = get(gca, 'YLim');
                                            line([13 13], [yBnds(1) yBnds(2)], 'Color', [.5 .5 .5])
                                            line([24 24], [yBnds(1) yBnds(2)], 'Color', [.5 .5 .5])
                                            ylabel('Unit ID');
                                            xlim([1 size(ppm,2)]);
                                            set(gca, 'XTick', 1:3:36, 'XTickLabel', {'T', 'F', 'P', 'R'});
                                            fix_font;
                                            
                                            
                                            
                                            
                                            subplot(3,1,3);
                                            hold on;
                                            imagesc(1:size(ppm,2), group(g).rat(r).day(d).binCtrs, ppm);
                                            line([allTms(1) allTms(end)], [180 180], 'LineWidth', 1.5, 'Color', [1 1 1]);
                                            axis xy;
                                            colormap jet;
                                            xlim([1 size(ppm,2)])
                                            ylim([0 360]);
                                            yBnds = get(gca, 'YLim');
                                            for c = 1:size(ppm,2)
                                                if ~isnan(maxInds(c))
                                                    plot(c, group(g).rat(r).day(d).binCtrs(maxInds(c)), '*w');
                                                end
                                            end
                                            plot(cycleInds, fitLineY, '--', 'Color', [.5 .5 .5]);
                                            line([13 13], [yBnds(1) yBnds(2)], 'Color', [.5 .5 .5])
                                            line([24 24], [yBnds(1) yBnds(2)], 'Color', [.5 .5 .5])
                                            ylabel('Position (degrees)');
                                            xlabel('Phase');
                                            fix_font;
                                            set(gca, 'XTick', 1:3:36, 'XTickLabel', {'T', 'F', 'P', 'R'})
                                            
                                        end %plotting each cycle
                                        
                                    end %if we can get enough data to decode
                                    
                                    
                                end %If # of cells in spkRstr with >=minNumSpks is >=minNumU, then decode
                                
                            end %runspeed & theta amplitude
                            
                            
                        end %cycles (peaks)
                        
                        
                    end
                    
                end %lap loop
            end %begin
            
        end %if there are more than 2 units
        
        
        %         figName = [group(g).name, '_Rat' num2str(r') '_Day', num2str(d)];
        figName = [group(g).name '_' group(g).rat(r).name '_' group(g).rat(r).day(d).name '_thetaTet' num2str(tetNum)];
        figure('name', figName);
        hold on;
        imagesc(1:36, group(g).rat(r).day(d).binCtrs, mean(dayPpm,3));
        line([1 36], [180 180], 'LineWidth', 1.5, 'Color', [1 1 1]);
        axis xy;
        colormap jet;
        xlim([1 36])
        ylim([0 360]);
        line([13 13], [0 360], 'Color', [.5 .5 .5])
        line([24 24], [0 360], 'Color', [.5 .5 .5])
        ylabel('Position (degrees)');
        set(gca, 'XTick', 1:3:36, 'XTickLabel', {'T', 'R', 'P', 'F'})
        xlabel('Phase');
        %         title([num2str(size(dayPpm,3)) ' Theta Cycles']);
        title({[group(g).rat(r).name]; [num2str(size(dayPpm,3)) ' Theta Cycles']; ['theta tet = csc' num2str(tetNum)]})
        fix_font;
        ylim([90 270]);
        colorbar;
        caxis([0.003 0.060001]);
        
        if saveFigs == 1
            curDir = pwd;
            cd(saveDir);
            %             savefig(figName)
            %             print(figName, '-dpng');
            saveas(gcf, figName, 'epsc')
            saveas(gcf, figName, 'png')
            saveas(gcf, figName, 'fig')
            cd(curDir);
        end
        
        group(g).rat(r).day(d).decErrs = decErr;
        group(g).rat(r).day(d).slopes = slopes;
        
    end %day
    
    %     end %rat
    
    
end %group



figure('Position', [ 420         362        1008         420]);
% slopeFigName = ['error_slope_tet' num2str(forFigName) '_day' num2str(3)];
for g = 1:2
    if g == 1
        r = 1;
        d = 1;
    else
        r = 2;
        d = 3;
    end
    %     subplot(2,2,(g-1)*2+1)
    subplot(1,2,1)
    hold on;
    %     histogram(group(g).rat(r).day(d).decErrs, 'BinWidth', 10, 'BinLimits', [0 180])
    try
        [cumProp, xScale] = calc_cum_prop_ip_range(group(g).rat(r).day(d).decErrs, [0 180]);
    catch
        keyboard
    end
    errFig(g) = plot(xScale, cumProp);%#ok
    if g == 2
        ylabel('Cumulative Proportion');
        xlabel('Decoding St Dev (degrees)');
        title('Decoding Error');
        fix_font;
        legend(errFig, {'WT', 'KO'}, 'Location', 'SouthEast');
        axis square;
    end
    
    %     subplot(2,2,(g-1)*2+2)
    subplot(1,2,2)
    hold on;
    %     histogram(group(g).rat(r).day(d).slopes, 'BinWidth', 5, 'BinLimits', [-50 50]);
    [cumProp, xScale] = calc_cum_prop_ip_range(group(g).rat(r).day(d).slopes, [-50 50]);
    slopeFig(g) = plot(xScale, cumProp); %#ok
    if g == 2
        ylabel('Cumulative Proportion');
        xlabel('Slope');
        title('Slope Distribution');
        fix_font;
        legend(slopeFig, {'WT', 'KO'}, 'Location', 'SouthEast');
        axis square;
    end
end


if saveFigs == 1
    curDir = pwd;
    cd(saveDir);
    %     savefig('DecodingErr_and_Slope_Distributions')
    %     print('DecodingErr_and_Slope_Distributions', '-dpng');
    saveas(gcf, slopeFigName, 'eps')
    saveas(gcf, slopeFigName, 'png')
    saveas(gcf, slopeFigName, 'fig')
    cd(curDir);
end


end %fnctn





