function fmr1CircTrack_x_thetaDecodeByTime(group)
% function fmr1CircTrack_x_thetaDecodeByTime(group)
%
% JB Trimper 01/2021 Colgin Lab

plotEachCycle = 1; %set to 1 to plot each cycle (waves, spikes, and PPM)
%               ****   If you choose to plot each cycle, you'll need to
%         enter a manual stop because there are hundreds and your PC will crash  ****

figCntr = 0; %initialize
maxNumFigs = 20; %stop point

saveFigs = 0; %Set to 1 to save figures
saveDir =  'E:\FMR1_CIRCTRACK\RESULTS\thetaCycleDecodingByTime';

shiftPlots = 1; %set to 1 to have individual cycle plots shifted so rat's current location is always plotted at 180
%                        ***This needs to be set to 1 for the average plot to be interpretable ***
%              *** Only set it to zero for the sake of visual inspection of individual theta cycle plots ***


% Decoding parameters
sampRate = 20000; %HzsampRate = 20000; %Hz - spike sampling rate
bayesWin = 0.0104; %s - decoder time window
bayesStep = 0.0104; %s - decoder step

minNumU = 3; %minimum # of cells that need to fire at least minNumSpks each
minNumSpks = 3; %these numbers from Zheng et al

minThetaAmp = -inf; %Minimum theta amplitude --> zscore((abs(hilbert(thetaFilteredLfp)))

runThresh = 5; %cm/s

decWinSamps = 750; %LFP samples - window around the middle of the theta cycle to do decoding
%                   = 0.375s (~3 cycles)

radBinCtrs = group(2).rat(1).day(1).binCtrs; %doesn't change across days/rats
newRewLoc = 0; %used to get all plots to line up as though rewards are at 0 & 180
[~,newRewInd] = min(abs(circ_dist(deg2rad(radBinCtrs), deg2rad(newRewLoc))-0));

curDir = pwd; %Note where we are so we can return to it later

meanAmps = [];

for g = 1:2
% for g = 1
    fprintf('Group %d\n', g);
    
    if g == 1
        r = 1;
    else
        r = 2;
    end
    
    fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);

    dayNums = 1:length(group(g).rat(r).day);
%     if g == 2
%         dayNums = 1; %So we're only considering the KO rat's day w/ most cells
%     else
%         dayNums = 1:2; 
%     end
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
        uIDs = zeros(length(group(g).rat(r).day(d).xBeginUnitInfo),2); %unit is bad if max firing rate in bin does not exceed 1
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
                
                lfpRoot = ['CSC' num2str(tetNum)];
                lfpStruct = read_in_lfp(['CSC' num2str(tetNum) '.ncs']);
                
                load([lfpRoot '_broadThetaLfp.mat']); %#ok
                lfpStruct.broadThetaLfp = filtLfp;
                
                load([lfpRoot '_narrowThetaLfp.mat']);  %#ok
                lfpStruct.narrowThetaLfp = filtLfp;
                thetaFiltLfpAmp = zscore(abs(hilbert(filtLfp)));
                
                cd(curDir);
                
                phiVctr = get_asym_theta_phi_vector(lfpStruct); %theta phase for each LFP sample
                phiTms = get_theta_phase_times(lfpStruct);
                
                %Identify rising phase to split theta time-series into cycles
                phiVctr = phiVctr+pi/2;
                zeroCrossBnry = diff(sign(phiVctr));
                risePhiInds = find(zeroCrossBnry>=1);
                risePhiInds(diff(risePhiInds)==1) = []; %get rid of immediate neighbors
                
                thetaCycleTimes = cell(1,length(risePhiInds)-1);
                for tc = 1:length(risePhiInds)-1
                    thetaCycleTimes{tc} = lfpStruct.ts(risePhiInds(tc):risePhiInds(tc+1));
                end
                
                
                lapInds = group(g).rat(r).day(d).begin(b).lapInds;
                lapTms = group(g).rat(r).day(d).begin(b).lapTms;
                
                
                
                for lp = 1:size(lapInds,1)
                    
                    startLpTm = lapTms(lp,1); %Referenced to video
                    endLpTm = lapTms(lp,2);
                    
                    % DEFINE EDGE OF CYCLE AS TROUGHS
                    lapPhiInds = find(phiTms(3,:)>=startLpTm & phiTms(3,:)<=endLpTm); %Use troughs -- find all troughs that occur during this lap
                    lapPhiTms = phiTms(3,lapPhiInds); %#ok
                    allInds = [allInds 1:length(lapPhiTms)-1];
                    
                    for p = 1:length(lapPhiTms)-1
                       
                        startCycle = lapPhiTms(p);
                        endCycle = lapPhiTms(p+1);
                        
                        midCycle = mean([startCycle endCycle]);
                        winStartTm = midCycle - ((decWinSamps/2)/2000);
                        winEndTm = midCycle + ((decWinSamps/2)/2000);
                        winDur = winEndTm - winStartTm;
                        
                        lfpInds = find(lfpStruct.ts>=startCycle & lfpStruct.ts<=endCycle);
                        meanAmp = mean(thetaFiltLfpAmp(lfpInds)); %#ok
                        meanAmps = [meanAmps meanAmp]; %#ok
                        
                        meanRs = mean(smRs(smRs(:,1)>=startCycle & smRs(:,1)<=endCycle,2));
                        if meanRs >= runThresh  &&  meanAmp >= minThetaAmp
                            
                            goodRs = [goodRs p];
                            
%                             if p == 2; keyboard; end
                            
                            %Pre-allocate for spike raster
                            nEvBins = round(winDur*sampRate);
                            spkRstr = zeros(size(uIDs,1), nEvBins);
                            
                            cycleSpkTms = cell(1,size(uIDs,1));
                            winSpkTms = cell(1,size(uIDs,1));
                            
                            uCntr = 0;
                            for u = 1:length(group(g).rat(r).day(d).begin(b).unit)
                                uID = group(g).rat(r).day(d).begin(b).unit(u).ID;
                                if ismember(uID, uIDs, 'row') %if the unit wasn't discarded due to low firing ratemap
                                    uCntr = uCntr + 1;
                                    
                                    allSpkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                                    
                                    tmpCycleSpkTms = allSpkTms(allSpkTms>=startCycle & allSpkTms<=endCycle);
                                    
                                    if length(tmpCycleSpkTms)>=minNumSpks %at least 3 spks
                                        
                                        goodSpks = [goodSpks p];
                                        
                                        cycleSpkTms{uCntr} = tmpCycleSpkTms;
                                        winSpkTms{uCntr} = allSpkTms(allSpkTms>=winStartTm & allSpkTms<=winEndTm);
                                        
                                        % Take spike times and fill in the raster
                                        timePassed = winSpkTms{uCntr} - winStartTm;
                                        spkInds = round(timePassed * sampRate);
                                        spkInds(spkInds==0)=1;
                                        
                                        spkRstr(uCntr, spkInds) = 1;
                                        
                                    end
                                    
                                end
                                
                            end %unit
                            
                            
                            
                            if sum(sum(spkRstr,2)>=minNumSpks)>=minNumU %If # of cells in spkRstr with >=minNumSpks is >=minNumU, then decode
                               
                                
                                %Plot the cycle
                                lfpInds = find(lfpStruct.ts>=winStartTm & lfpStruct.ts<=winEndTm);
                                if length(lfpInds)>750
                                    lfpInds(751:end) = [];
                                end
                                if lfpInds(1)>0   &&   lfpInds(end)<length(lfpStruct.ts)
                                    
                                goodUnits = [goodUnits p];
                                
                                
                                    goodLfpInds = [goodLfpInds p];
                                    
                                    % Do the decoding
                                    ppm = BayesianDecoder(spkRstr,rateMaps,bayesWin,bayesStep,sampRate);
                                    ppm(isnan(ppm)) = 1/size(rateMaps,2);
                                    posInds = find(radPos(:,1)>=winStartTm & radPos(:,1)<=winEndTm);
                                    posInds = [posInds(1)-1; posInds; posInds(end)+1]; %#ok
                                    tmInds = 1:size(ppm,2);
                                    winMids = winStartTm+((tmInds-1)*bayesStep)+(.5*bayesWin);
                                    
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
                                    cycleInds = find(winMids>=startCycle & winMids<=endCycle);
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
                                        dayLfpData = lfpStruct.data(lfpInds);
                                    else
                                        dayPpm = cat(3,dayPpm, ppm);
                                        dayLfpData = [dayLfpData lfpStruct.data(lfpInds(1):lfpInds(end))]; %#ok
                                    end
                                    
                                    
                                    
                                    if plotEachCycle == 1 && figCntr < maxNumFigs
                                        keyboard
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
                                        
                                        
                                        subplot(3,1,2);
                                        hold on;
                                        uCntr = 0;
                                        for u = 1:length(cycleSpkTms)
                                            if ~isempty(cycleSpkTms{u})
                                                uCntr = uCntr + 1;
                                                for st = 1:length(winSpkTms{u})
                                                    line([winSpkTms{u}(st) winSpkTms{u}(st)], [uCntr-.3 uCntr+.3])
                                                end
                                            end
                                        end
                                        ylim([0 uCntr+1]);
                                        yBnds = get(gca, 'YLim');
                                        line([startCycle startCycle], [yBnds(1) yBnds(2)], 'Color', [.5 .5 .5])
                                        line([endCycle endCycle], [yBnds(1) yBnds(2)], 'Color', [.5 .5 .5])
                                        xlim(lfpStruct.ts([lfpInds(1) lfpInds(end)]))
                                        ylabel('Unit ID');
                                        fix_font;
                                        
                                        
                                        subplot(3,1,3);
                                        hold on;
                                        imagesc(winMids, group(g).rat(r).day(d).binCtrs, ppm);
                                        if shiftPlots == 0
                                            plot(winMids, allPos, 'LineWidth', 1.5, 'Color', [1 1 1]);
                                        else
                                            line([allTms(1) allTms(end)], [180 180], 'LineWidth', 1.5, 'Color', [1 1 1]);
                                        end
                                        axis xy;
                                        colormap jet;
                                        xlim(lfpStruct.ts([lfpInds(1) lfpInds(end)]))
                                        ylim([0 360]);
                                        yBnds = get(gca, 'YLim');
                                        for c = 1:size(ppm,2)
                                            if ~isnan(maxInds(c))
                                                plot(winMids(c), group(g).rat(r).day(d).binCtrs(maxInds(c)), '*w');
                                            end
                                        end
                                        plot(winMids(cycleInds), fitLineY, '--', 'Color', [.5 .5 .5]);
                                        line([startCycle startCycle], [yBnds(1) yBnds(2)], 'Color', [.5 .5 .5])
                                        line([endCycle endCycle], [yBnds(1) yBnds(2)], 'Color', [.5 .5 .5])
                                        ylabel('Position (degrees)');
                                        xlabel('Time (s)');
                                        fix_font;
                                    end %plotting each cycle
                                   
                                end %if we can get enough data to decode
                                
                                
                                    
                            end %If # of cells in spkRstr with >=minNumSpks is >=minNumU, then decode
                           
                        end %runspeed & theta amplitude
                        
                        
                    end %cycles (peaks)
                    
                end %lap loop
            end %begin
            
            
%             keyboard
            
        end %if there are more than 2 units
        
%         figName = [group(g).name, '_Rat' num2str(r') '_Day', num2str(d)];
        figName = [group(g).name '_' group(g).rat(r).name '_' group(g).rat(r).day(d).name];
        figure('name', figName);
        hold on;
        ppmTmVals = winMids - winMids(median([1:size(winMids,2)])); %#ok
        imagesc(ppmTmVals, group(g).rat(r).day(d).binCtrs, mean(dayPpm,3));
        line([ppmTmVals(1) ppmTmVals(end)], [180 180], 'LineWidth', 1.5, 'Color', [1 1 1]);
        axis xy;
        colormap jet;
        xlim([ppmTmVals(1) ppmTmVals(end)])
        ylim([0 360]);
        lfpTmVals = linspace(ppmTmVals(1), ppmTmVals(end), size(dayLfpData,1));
        avgLfp = mean(dayLfpData,2);
        avgLfp = rescale(avgLfp, 100, 130);
        plot(lfpTmVals, avgLfp, 'Color', [.8 .8 .8]);
        ylabel('Position (degrees)');
        xlabel('Time (s)');
        title([num2str(size(dayPpm,3)) ' Theta Cycles']);
        fix_font;
        colorbar;
        ylim([90 270]);
%         caxis([.007 0.042])
        
        if saveFigs == 1
            curDir = pwd;
            cd(saveDir);
            savefig(figName)
            print(figName, '-dpng');
            cd(curDir);
        end
        
        group(g).rat(r).day(d).decErrs = decErr;
        group(g).rat(r).day(d).slopes = slopes;
        
    end %day
    
    %     end %rat
    
    
end %group
keyboard


figure('Position', [ 420         362        1008         420]);
for g = 1:2
    if g == 1
        r = 1;
        d = 2;
    else
        r = 2;
        d = 2;
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
    savefig('DecodingErr_and_Slope_Distributions')
    print('DecodingErr_and_Slope_Distributions', '-dpng');
    cd(curDir);
end



end %fnctn





