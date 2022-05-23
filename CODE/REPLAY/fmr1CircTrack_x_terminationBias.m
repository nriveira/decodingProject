function fmr1CircTrack_x_terminationBias(group)
% function fmr1CircTrack_x_terminationBias(group)
%
%   PURPOSE:
%       For comparing the termination and initiation bias at reward
%       location for WT and KO groups.
%
%   INPUTS:
%       group struct, post-function 6.
%
%   OUTPUTS:
%
%   OPTIONS:
%       r2Thresh = threshold for r2 values for event to be considered a
%           replay event
%       saveOrNot (0 = don't save, 1 = save figs)
%
% MMD
% Colgin Lab
% 06/2021

%% OPTIONS

r2Thresh = 0.50; %as in Zheng et al. 2021

plotEachDay = 0;
saveDayPlots = 0;

saveOrNot = 1; %to save figs that aren't day plots

saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\REPLAY\terminationBias';
spatBinSz = 4;
% cols = {'Blue', 'Red'};
% nboot = 5000; %for calculating confidence intervals

%% INITIALIZE

% Decoding parameters
sampRate = 20000; %HzsampRate = 20000; %Hz - spike sampling rate
bayesWin = 20/1000; %20 ms time window, as in Hwaun & Colgin 2019
bayesStep = 10/1000; %10 ms time step, as in Hwaun & Colgin 2019

numTimeBins = 100;
timeBinVals = 0:1/numTimeBins:1;
numPosBins = 360/spatBinSz;

ctrlPosxTime = cell(2,2); %by group, before learning reward location for the day (aka sleep1)
posxTime = cell(2,1); %for storing by group
posxTimeBySlope = cell(2,2); %group x slope (forward or reverse event)

radBinCtrs = group(2).rat(1).day(1).binCtrs; %doesn't change across days/rats
radBinCtrs = deg2rad(radBinCtrs);

newRewLoc = 0; %used to get all plots to line up as though rewards are at 0 & 180
[~,newRewInd] = min(abs(circ_dist(deg2rad(radBinCtrs), deg2rad(newRewLoc))-0));

curDir = pwd; %for returning after saving figs
if saveOrNot == 1 || saveDayPlots == 1
    cd(saveDir)
end

ttls = {'All events', 'Forward events', 'Reverse events'};

%% GET DATA

for g = 1:2
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            if ~isempty(group(g).rat(r).day(d).sleep(2).coords) %if there is any sleep data
                dayPosxTime = cell(1,3); %all, for, rev
                
                rewLocs = group(g).rat(r).day(d).rewLocs;
                [~,rewInd] = min(abs(circ_dist(deg2rad(radBinCtrs), deg2rad(rewLocs(1)))-90)); %shift so reward locations are always 90 and 270
                shiftVal = newRewInd - rewInd;
                
                rateMaps = zeros(length(group(g).rat(r).day(d).xAllBeginUnitInfo), length(group(g).rat(r).day(d).xAllBeginUnitInfo(1).smRateMap));
                
                badU = [];
                uIDs = zeros(length(group(g).rat(r).day(d).xAllBeginUnitInfo),2);
                for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                    if max(group(g).rat(r).day(d).xAllBeginUnitInfo(u).rateMap)>=1 %unit is bad if max firing rate in bin does not exceed 1
                        rateMaps(u,:) = group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap; %Smoothed ratemap
                        uIDs(u,:) = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID;
                    else
                        badU = [badU u]; %#ok
                    end
                end
                rateMaps(badU,:) = [];
                rateMaps(rateMaps==0) = 0.0001; %get rid of zeros because our Bayesian decoder can't handle 'em.
                uIDs(badU,:) = [];
                
                %             rateMaps = circshift(rateMaps, [0 shiftVal]); going to circ
                %             shift ppm later
                
                %             for s = 1:5 %Need to SPLIT off s1 later!!!!
                for s = 2:5
                    if ~isempty(group(g).rat(r).day(d).sleep(s).unit)
                        popEvents = group(g).rat(r).day(d).sleep(s).popEvents; %shorten variable name
                        
                        for i = 1:length(popEvents)
                            startTm = popEvents(i,1);
                            endTm = popEvents(i,2);
                            
                            nEvBins = round((endTm-startTm) * sampRate); %number of bins in spike raster
                            spkRstr = zeros(size(uIDs,1), nEvBins);
                            
                            timeaxis = 0:bayesStep:(endTm-startTm);
                            
                            uCntr = 0; %can't just use u since we discard some units due to firing rate
                            for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                                
                                uID = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID;
                                if ismember(uID, uIDs, 'row') %if the unit wasn't discarded due to low firing ratemap
                                    uCntr = uCntr + 1;
                                    
                                    allSpkTms = group(g).rat(r).day(d).sleep(s).unit(u).spkTms;
                                    tmpEvSpkTms = allSpkTms(allSpkTms>=startTm & allSpkTms<=endTm);
                                    
                                    % Take spike times and fill in the raster
                                    timePassed = tmpEvSpkTms - startTm;
                                    spkInds = round(timePassed * sampRate);
                                    spkInds(spkInds==0)=1;
                                    
                                    spkRstr(uCntr, spkInds) = 1;
                                    
                                end %if unit wasn't discarded
                            end %units
                            
                            ppm = BayesianDecoder(spkRstr,rateMaps,bayesWin,bayesStep,sampRate); %Ernie's decoder
                            ppm(isnan(ppm)) = 1/size(rateMaps,2); %nan where pxn is chance
                            
                            [maxVals, maxInds] = max(ppm);
                            maxInds(maxVals==1/size(rateMaps,2)) = NaN;
                            
                            bins2use = 1:size(ppm,2);
                            bins2use(isnan(maxInds)) = []; %don't use these
                            
                            if length(bins2use)/size(ppm,2) > .8 %as long as spikes in at least 80% of bins - I don't think find_population_events does this as of writing this code but it should
                                [r2, ~, ~, ~, slope] = Cir_reg(ppm, radBinCtrs, timeaxis, bins2use); %CZ code
                                
                                if slope > 0
                                    sInd = 1; %positive slope
                                else
                                    sInd = 2;
                                end %forward or reverse event
                            else
                                r2 = nan;
                            end
                            
                            if r2 > r2Thresh
                                %need to normalize time to 0-1
                                
                                tmpPosxTime = zeros(numPosBins,numTimeBins);
                                
                                evDur = endTm - startTm;
                                ppmTAxis = 0:evDur/ size(ppm,2):evDur;
                                tmpBinTm = rescale(ppmTAxis, 0, 1);
                                
                                newBinStart = 1;
                                for bn = 2:length(tmpBinTm)
                                    newBinEnd = round(tmpBinTm(bn),2);
                                    newBinEnd = round(newBinEnd/(1/numTimeBins),0);
                                    try
                                        tmpPosxTime(:,newBinStart:newBinEnd) = repmat(ppm(:,bn-1), [1, newBinEnd-newBinStart+1]); %tile to fit norm time distribution
                                    catch; keyboard; end
                                    newBinStart = newBinEnd;
                                end %bin
                                
                                dayPosxTime{1} = cat(3, dayPosxTime{1}, tmpPosxTime);
                                dayPosxTime{sInd+1} = cat(3, dayPosxTime{sInd+1}, tmpPosxTime);
                                
                                shiftPosxTime = circshift(tmpPosxTime, shiftVal, 1);
                                
                                if s == 1
                                    ctrlPosxTime{g} = cat(3, ctrlPosxTime{g}, shiftPosxTime);
                                else
                                    posxTime{g} = cat(3, posxTime{g}, shiftPosxTime);
                                    posxTimeBySlope{g,sInd} = cat(3, posxTimeBySlope{g,sInd}, shiftPosxTime);
                                end
                                
                                
                            end %r2 is higher than threshold set
                        end %pop events
                        
                    end %if there are pop events
                end %sleep
                
                if plotEachDay == 1
                    figtitle = ['terminationBias_' group(g).rat(r).name '_' group(g).rat(r).day(d).name];
                    
                    figure('Name', figtitle, 'Position', [477 596 1048 385])
                    
                    tmpMax = 0;
                    for sp = 1:3
                        pullPxns = dayPosxTime{sp};
                        meanPxn = mean(pullPxns,3);
                        if max(meanPxn(:)) > tmpMax
                            tmpMax = max(meanPxn(:));
                        end
                        subplot(1,3,sp)
                        imagesc(timeBinVals, rad2deg(radBinCtrs), meanPxn)
                        ylim([0 360])
                        axis xy
                        axis square
                        colormap parula
                        
                        for rw = 1:length(group(g).rat(r).day(d).rewLocs)
                            if group(g).rat(r).day(d).rewLocs(rw) ~= 0
                                line([0 1], [group(g).rat(r).day(d).rewLocs(rw) group(g).rat(r).day(d).rewLocs(rw)], 'Color', 'white', 'LineStyle', '--', 'LineWidth', 1.5)
                            else
                                line([0 1], [2 2], 'Color', 'white', 'LineStyle', '--', 'LineWidth', 1.5)
                                line([0 1], [358 358], 'Color', 'white', 'LineStyle', '--', 'LineWidth', 1.5)
                            end
                        end %reward
                        
                        line([0 1], [90 90], 'Color', 'black', 'LineStyle', '--', 'LineWidth', 1.5)
                        xticks(0:0.5:1)
                        yticks(0:90:360)
                        
                        title(ttls{sp})
                    end %subplot
                    
                    caxis([0 tmpMax])
                    cbr = colorbar;
                    cbr.Position = [0.92 0.3169 0.0220 0.4000];
                    ylabel(cbr, 'Probability')
                    
                    if saveDayPlots == 1
                        saveas(gcf, figtitle, 'epsc');
                        saveas(gcf, figtitle, 'fig');
                        saveas(gcf, figtitle, 'png');
                    end %save day plots
                end %plot each day
            end %if there is any sleep data
        end %day
        
    end %rat
end %group

keyboard

%% FIG 1

figtitle = 'TerminationBias';

figure('Name', figtitle, 'Position', [428 283 912 707])
meanForPlot = zeros(90,100,2,3); %pos x time x group x slope (all, forward, reverse)

for g = 1:2
    pullMaps = posxTime{g};
    meanForPlot(:,:,g,1) = mean(pullMaps, 3);
    
    for sInd = 1:2
        pullMaps = posxTimeBySlope{g,sInd};
        meanForPlot(:,:,g,sInd+1) = mean(pullMaps,3);
    end %slope ind
end %group

maxAll = max(meanForPlot(:)); %so same color scale can be used for both
spMap = [1 2; 3 4; 5 6];
for g = 1:2
    for row = 1:3 %row of subplot
    subplot(3,2,spMap(row,g))
    
    imagesc(timeBinVals, rad2deg(radBinCtrs), meanForPlot(:,:,g,row))
    axis xy
    axis square
    
    colormap(parula)
    caxis([0 maxAll])
    freezeColors
    
    line([0 100], [90 90], 'Color', 'white', 'LineStyle', '--', 'LineWidth', 1.5) %make line at 90
    line([0 100], [270 270], 'Color', 'white', 'LineStyle', '--', 'LineWidth', 1.5) %make line at 270
    
    xlim([0 1])
    xticks([0 0.5 1])
    if row == 3
    xlabel('Normalized time in SWR')
    end
    
    ylim([1 360])
    yticks(0:90:360)
    ylabel('Angular position (deg)')
    
    title([group(g).name ' - ' ttls{row}])
    end %subplot row
end %group

cbr = colorbar;
cbr.Position = [0.85 0.71 0.0248 0.2136];
ylabel(cbr, 'Probability')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %save day plots

cd(curDir)

end %function


