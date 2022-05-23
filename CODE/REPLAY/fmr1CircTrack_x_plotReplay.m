function fmr1CircTrack_x_plotReplay(group)
% function fmr1CircTrack_x_plotReplay(group)
%
%   PURPOSE:
%       For comparing the replay qualities from the sleeps between WT and
%       KO groups.
%
%   INPUTS:
%       group struct, post-function 6.
%
%   OUTPUTS:
%       F1: Distribution of all r2 values.
%       F2: Distribution of r2 values, separated by forward and reverse
%           events.
%       F3: Pie charts of forward vs reverse events (all and R^2 > 0.5).
%       F4: Absolute value of the slopes of all replay events (R^2 > 0.5).
%       F5: Absolute value of the slopes of replay events (R^2 > 0.5),
%           separated by forward and reverse events.
%       F6: Plot replay event duration (all events).
%       F7: Plot and calculated correlationof event duration x normlized
%           slope.
%       F8: Plot of 5 highest r2 values from both groups (optional).
%       F9: Plot of 5 random r2 values that fall between a set min and max
%           value for both groups (optional) - just in case you want more
%           examples. Can set other qualifications too
%
%   OPTIONS:
%       saveOrNot: save figs (1) or don't (0)
%       See function for other options
%
% MMD
% Colgin Lab
% 06/2021

%% OPTIONS

saveOrNot = 1; %to save figs

prepForStats = 0; %1 to prepare output for SPSS and keyboard in code

downSampEvents = 1; %to randomly down sample events to to same as lowest group
if downSampEvents == 1
    rsE = RandStream('mt19937ar');
end

downSampCell = 0; %1 to downsamp, 0 to not
if downSampCell == 1
    newCellNum = 40; %number of cells to down sample to
    rsC = RandStream('mt19937ar');
end

minFr = 1; %Hz, for cell to be included

plotHighest = 1; %whether plot the 5 highest r2 values for each group (1 = plot)

plotRandom = 0; %whether plot first 5 replay sequences (r2>0.7 or whatever below) for each group (1 = plot)
%  Check lines 231-236 for how it actually makes the selection.
randMin = 0.6; %min for random r2 to plot
randMax = 1; %change to 1 if you don't want a max
randSumMaxVals = .6; %sum up for that higher probability decoded event

saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\REPLAY\replayProperties';

cols = {'Blue', 'Red'};
nboot = 5000; %for calculating confidence intervals

jitter = 0.2; %for dotplots

r2Thresh = 0.5;

%% INITIALIZE

% Decoding parameters
sampRate = 20000; %HzsampRate = 20000; %Hz - spike sampling rate
bayesWin = 20/1000; %20 ms time window, as in Hwaun & Colgin 2019
bayesStep = 10/1000; %10 ms time step, as in Hwaun & Colgin 2019

rhoValues = cell(2,1); %for storing r2 values
rhoValuesBySlope = cell(2,2); %group x forward or reverse

ForRevAll = cell(2,1); %this seems like a weird way to do it but I need to to downsample
ForRevR2 = cell(2,1); %for events that meet R2 > 0.5 threshold

slopesAll = cell(2,1); %absolute values of all slopes
slopesbyEvType = cell(2,2); %forwards vs. reverse replay events

normSlopesAll = cell(2,1); %rescaled so time of SWR event is 0-1
normSlopesbyEvType = cell(2,2);

avgJumpAll = cell(2,1); % average spatial jump for decoded positions in replay event
avgJumpbyEvType = cell(2,2);

pathDistAll = cell(2,1); %just dist between end and beginning of decoded path
pathDistbyEvType = cell(2,2);

eventDurs = cell(2,1); %for storing event durations
evDursbyEvType = cell(2,2); %forwards vs. reverse

evDurxDist = cell(2,1); %for event duration x trajectory length

radBinCtrs = group(2).rat(1).day(1).binCtrs; %doesn't change across days/rats
radBinCtrs = deg2rad(radBinCtrs);

slopeNames = {'Forward', 'Reverse'};

if plotHighest == 1
    ppmForPlot = cell(5,2); %for storing the ppms for plotting
    r2ForPlot = zeros(5,2); %for storing the matching r2 values
end

if plotRandom == 1
    ppmForPlotRand = cell(5,2); %for storing the ppms for plotting
    r2ForPlotRand = zeros(5,2); %for storing the matching r2 values
    meanProb = zeros(5,2);
end %plot random

curDir = pwd; %for returning after saving figs

%% GET DATA

for g = 1:2
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            rateMaps = zeros(length(group(g).rat(r).day(d).xAllBeginUnitInfo), length(group(g).rat(r).day(d).xAllBeginUnitInfo(1).smRateMap));
            
            badU = [];
            uIDs = zeros(length(group(g).rat(r).day(d).xAllBeginUnitInfo),2);
            for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                if max(group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap) >= minFr
                    rateMaps(u,:) = group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap; %Smoothed ratemap
                    uIDs(u,:) = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID;
                else
                    badU = [badU u]; %#ok
                end %if meets max
            end %units
            rateMaps(badU,:) = [];
            rateMaps(rateMaps==0) = 0.0001; %get rid of zeros because our Bayesian decoder can't handle 'em.
            uIDs(badU,:) = [];
            
            
            if downSampCell == 1 && size(rateMaps,1) > newCellNum
                y = datasample(rsC,1:size(rateMaps,1),newCellNum,'Replace',false);
                y = sort(y);
                
                uIDs = uIDs(y,:);
                rateMaps = rateMaps(y,:);
            end %down sample cells
            
            for s = 2:length(group(g).rat(r).day(d).sleep) %after experience, so just sleep 2 on
                if ~isempty(group(g).rat(r).day(d).sleep(s).unit) %if there is anything in the sleep
                    
                    popEvents = group(g).rat(r).day(d).sleep(s).popEv; %shorten variable name
                    
                    for i = 1:length(popEvents)
                        startTm = popEvents(i).tms(1);
                        endTm = popEvents(i).tms(2);
                        
                        eventDurs{g} = [eventDurs{g} endTm-startTm];
                        
                        r2 = popEvents(i).r2;
                        slope = popEvents(i).slope;
                        pxn = popEvents(i).pxn;
                        
                         tmpForRev = zeros(1,2);
                        if slope > 0 %forward event
                            sInd = 1; %positive slope
                            tmpForRev(1) = 1;
                            
                        else %negative event
                            sInd = 2;
                            tmpForRev(2) = 1;
                        end %which direction is slope
                        
%                         nEvBins = round((endTm-startTm) * sampRate); %number of bins in spike raster
%                         spkRstr = zeros(size(uIDs,1), nEvBins);
%                         
%                         timeaxis = 0:bayesStep:(endTm-startTm);
%                         
%                         uCntr = 0; %can't just use u since we discard some units due to firing rate (or down sampling)
%                         for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
%                             
%                             uID = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID;
%                             if ismember(uID, uIDs, 'row') %if the unit wasn't discarded due to low firing ratemap
%                                 uCntr = uCntr + 1;
%                                 
%                                 allSpkTms = group(g).rat(r).day(d).sleep(s).unit(u).spkTms;
%                                 tmpEvSpkTms = allSpkTms(allSpkTms>=startTm & allSpkTms<=endTm);
%                                 
%                                 % Take spike times and fill in the raster
%                                 timePassed = tmpEvSpkTms - startTm;
%                                 spkInds = round(timePassed * sampRate);
%                                 spkInds(spkInds==0)=1;
%                                 
%                                 spkRstr(uCntr, spkInds) = 1;
%                                 
%                             end %if unit wasn't discarded
%                         end %units
                        
%                         ppm = BayesianDecoder(spkRstr,rateMaps,bayesWin,bayesStep,sampRate); %Ernie's decoder
                        pxn(isnan(pxn)) = 1/size(rateMaps,2); %nan where pxn is chance - change for plotting purposes
                        
                        [maxVals, maxInds] = max(pxn);
                        maxInds(maxVals==1/size(rateMaps,2)) = NaN;
%                         
%                         bins2use = 1:size(ppm,2);
%                         bins2use(isnan(maxInds)) = []; %don't use these - not better than chance
%                         
%                         if length(bins2use)/size(ppm,2) >= .8 %as long as spikes in at least 80% of bins - I don't think find_population_events does this as of writing this code but it should
%                             [r2,~,~,~,slope]  = Cir_reg(ppm, radBinCtrs, timeaxis, bins2use); %CZ code - slope in rad/s
%                             
%                             tmpForRev = zeros(1,2);
%                             
%                             
%                         else %more than 20% bins no spike
%                             r2 = nan;
%                         end %no spike bins
                        
                        if ~isnan(r2) %just get rid of nans here
                            ForRevAll{g} = [ForRevAll{g}; tmpForRev];
                            if r2 >= r2Thresh
                                ForRevR2{g} = [ForRevR2{g}; tmpForRev];
                                
                                evDursbyEvType{g,sInd} = [evDursbyEvType{g,sInd} endTm-startTm];
                            end %replay event (r2 over thresh)
                            
                            rhoValues{g} = [rhoValues{g} r2];
                            rhoValuesBySlope{g,sInd} = [rhoValuesBySlope{g,sInd} r2];
                            
                            if r2 >= r2Thresh
                                slopesAll{g} = [slopesAll{g} rad2deg(abs(slope))];
                                slopesbyEvType{g,sInd} = [slopesbyEvType{g,sInd} rad2deg(abs(slope))];
                                
                                %slope is in rad/s, need to convert to rad/normalized replay event
                                eventDur = endTm-startTm;
                                normSlope = slope .* eventDur; %rescale slope
                                
                                normSlopesAll{g} = [normSlopesAll{g} rad2deg(abs(normSlope))];
                                normSlopesbyEvType{g,sInd} = [normSlopesbyEvType{g,sInd} rad2deg(abs(normSlope))];
                                
                                notNan = find(~isnan(maxInds));
                                tmpDist = circ_dist(radBinCtrs(maxInds(notNan(end))), radBinCtrs(maxInds(notNan(1)))); %distance between two points
                                
                                if sInd == 1 && tmpDist < 0  %circ_dist just takes the shortest abs val amount, + or -
                                    tmpDist = 2*pi + tmpDist;
                                elseif sInd == 2 && tmpDist > 0
                                    tmpDist = 2*pi - tmpDist;
                                end
                                pathDistAll{g} = [pathDistAll{g} rad2deg(abs(tmpDist))];
                                pathDistbyEvType{g,sInd} = [pathDistbyEvType{g,sInd} rad2deg(abs(tmpDist))];
                                
                                jumpDist = nan(1,length(maxInds)-1);
                                for ji = 1:length(jumpDist)
                                    if ~isnan(maxInds(ji+1)) && ~isnan(maxInds(ji))
                                        tmpJump = circ_dist(radBinCtrs(maxInds(ji+1)), radBinCtrs(maxInds(ji)));
                                        jumpDist(ji) = abs(tmpJump);
                                    end
                                end %jumpInd
                                
                                avgJumpAll{g} = [avgJumpAll{g} rad2deg(nanmean(jumpDist))];
                                avgJumpbyEvType{g,sInd} = [avgJumpbyEvType{g,sInd} rad2deg(nanmean(jumpDist))];
                                
                                evDurxDist{g} = [evDurxDist{g}; eventDur rad2deg(abs(tmpDist))]; %also get event dur x normalized slope
                            end %r2 for slope
                        end %r2 is not nan
                        
                        if plotHighest == 1 && r2 > min(r2ForPlot(:,g))
                            [~, storeInd] = min(r2ForPlot(:,g));
                            ppmForPlot{storeInd,g} = pxn;
                            r2ForPlot(storeInd,g) = r2;
                        end %if we need to store this for plotting
                        
                        if plotRandom == 1 && r2 > randMin && r2 < randMax && r2 > min(r2ForPlotRand(:,g)) && sum(maxVals) > randSumMaxVals
                            % if plotRandom == 1 && mean(maxVals) > min(meanProb(:,g)) && r2 > randMin
                            [~, storeInd] = min(meanProb(:,g));
                            ppmForPlotRand{storeInd,g} = pxn;
                            r2ForPlotRand(storeInd,g) = r2;
                            
                            meanProb(storeInd,g) = sum(maxVals);
                            
                        end %plot random
                    end %pop events
                end %if there are pop events
            end %sleep
        end %day
    end %rat
end %group
keyboard


%% FIG 1 - REPLAY FIDELITY BY GROUP

figtitle = 'ReplayFidelity';

if downSampEvents == 1
    minEv = min([length(rhoValues{1}) length(rhoValues{2})]);
    
    for g = 1:2
        if length(rhoValues{1}) ~= minEv
            y = datasample(rsE,1:length(rhoValues{g}),minEv,'Replace',false);
            y = sort(y);
            rhoValues{g} = rhoValues{g}(y);
        end %if we have to downsample this one
    end %group
    
end %we have to down sample events

if downSampCell == 1
    figtitle = [figtitle '_downSampled_' num2str(newCellNum) 'cell'];
end
if downSampEvents == 1
    figtitle = [figtitle '_downSampled_events'];
end
figure('Position', [717 504 525 388], 'Name', figtitle)

lh = zeros(1,2);
leg = cell(1,2);

for g = 1:2
    rhoInGroup = rhoValues{g};
    
    nSWR = length(rhoInGroup);
    
    [~, sortInd] = sort(rhoInGroup, 'descend');
    r2Sorted = rhoInGroup(sortInd)';
    %     r2min = min(r2Sorted);
    r2min = 0;
    %     r2max = max(r2Sorted);
    r2max = 1;
    sigma = 0.1; %same as in Zheng et al. 2021
    
    [Wdis_pre,axisValue1] = WeightedProportion(r2Sorted, r2min, r2max, sigma); %from Ernie's code for Zheng et al. 2021
    
    sample = bootstrp(nboot, @WeightedProportion, r2Sorted, r2min, r2max, sigma);
    [CI_u, CI_l] = CorrectionForMultipleComparsion(sample);
    CIforPlot = [CI_u'; CI_l'];
    
    lh(g) = plot_filled_ci(axisValue1, Wdis_pre, CIforPlot, rgb(cols{g}));
    
    leg{g} = [group(g).name ' (n = ' num2str(nSWR) ')'];
    
end %group

ylabel('Proportion')
ylim([0 0.05])
yticks(0:0.01:0.05)

xlabel('r^2')

legend(lh, leg)

if saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    cd(curDir)
end

%% FIG 2 - REPLAY FIDELITY BY GROUP AND SLOPE

figtitle = 'ReplayFidelity_byForwardReverse';

if downSampCell == 1
    figtitle = [figtitle '_downSamp_' num2str(newCellNum) 'cell'];
end
if downSampEvents == 1
    figtitle = [figtitle '_downSampled_events'];
end
figure('Position', [553 487 910 388], 'Name', figtitle)

lh = zeros(2,2);
leg = cell(2,2);

for sInd = 1:2 %slope ind
    subplot(1,2,sInd)
    
    if downSampEvents == 1
        minEv = min([length(rhoValuesBySlope{1,sInd}) length(rhoValuesBySlope{2,sInd})]);
    end
    
    for g = 1:2
        
        rhoInGroupSlope = rhoValuesBySlope{g,sInd};
        if downSampEvents == 1 && length(rhoInGroupSlope) ~= minEv
            
            y = datasample(rsE,1:length(rhoInGroupSlope),minEv,'Replace',false);
            y = sort(y);
            rhoInGroupSlope = rhoInGroupSlope(y);
        end
        
        nSWR = length(rhoInGroupSlope);
        
        [~, sortInd] = sort(rhoInGroupSlope, 'descend');
        r2Sorted = rhoInGroupSlope(sortInd)';
        %         r2min = min(r2Sorted);
        %         r2max = max(r2Sorted);
        r2min = 0;
        r2max = 1;
        sigma = 0.1; %same as in Zheng et al. 2021
        
        [Wdis_pre,axisValue1] = WeightedProportion(r2Sorted, r2min, r2max, sigma); %from Ernie's code for Zheng et al. 2021
        
        sample = bootstrp(nboot, @WeightedProportion, r2Sorted, r2min, r2max, sigma);
        [CI_u, CI_l] = CorrectionForMultipleComparsion(sample);
        CIforPlot = [CI_u'; CI_l'];
        
        lh(g,sInd) = plot_filled_ci(axisValue1, Wdis_pre, CIforPlot, rgb(cols{g}));
        
        leg{g,sInd} = [group(g).name ' (n = ' num2str(nSWR) ')'];
        
    end %group
    
    ylabel('Proportion')
    xlabel('r^2')
    
    if isempty(leg{2}) %right now I have nothing for KO
        legend(lh(1,sInd), leg{1,sInd})
    else
        
        legend(lh(:,sInd), leg{:,sInd}) %fix this when I have KO data
    end
    
    %     ANOVA? for later
    %     [P,~,STATS] = ranksum(rhoValuesBySlope{1,sInd}, rhoValuesBySlope{2,sInd});
    %     fprintf('Replay fidelity - %s\n', slopeNames{sInd})
    %     fprintf('\tMann-Whitney U = %.2g, n1 = %d, n2 = %d, p = %.3g\n', STATS.ranksum, length(rhoValues{1,1}), length(rhoValues{2,1}), P);
    %
    %     if P < 0.05
    %         xlabel({'r^2'; ['Mann-Whitney U Test: p = ' num2str(P)]})
    %     end
    
    title(slopeNames{sInd})
    ylim([0 0.05])
    yticks(0:0.01:0.05)
    
end %sInd

if saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    cd(curDir)
end

%% FIG 3 - FORWARD VS REVERSE %

figtitle = 'Forward_reverse_Percentages';
if downSampCell == 1
    figtitle = [figtitle '_downSamp_' num2str(newCellNum) 'cell'];
end
if downSampEvents == 1
    figtitle = [figtitle '_downSampled_events'];
end

figure('Name', figtitle, 'Position', [419 457 1078 420])

labels = {'Forward', 'Reverse'};
pieCols = {'Purple', 'Green'};
spMap = [1 2; 3 4];
% ttlAdd = {' - All detected events', [' - Replay events (R^2 > ' num2str(r2thresh) ')']};

%first - all events

if downSampEvents == 1
    minEv = min([size(ForRevAll{1},1) size(ForRevAll{2},1)]);
    
    for g = 1:2
        if size(ForRevAll{g},1) ~= minEv
            y = datasample(rsE,1:size(ForRevAll{g},1),minEv,'Replace',false);
            y = sort(y);
            ForRevAll{g} = ForRevAll{g}(y,:);
        end %if this is the one we down sample
    end %group
    
    minEv = min([size(ForRevR2{1},1) size(ForRevR2{2},1)]);
    
    for g = 1:2
        if size(ForRevR2{g},1) ~= minEv
            y = datasample(rsE,1:size(ForRevR2{g},1),minEv,'Replace',false);
            y = sort(y);
            ForRevR2{g} = ForRevR2{g}(y,:);
        end %if this is the one we down sample
    end %group
    
end
rInd = 1;

for g = 1:2
    subplot(2,2,spMap(rInd,g))
    
    sumEv = sum(ForRevAll{g});
    if g == 1
        contData = sumEv;
    else
        koData = sumEv;
    end %for chi squ test
    h = pie([sumEv(1) sumEv(2)]);
    
    patchHand = findobj(h, 'Type', 'Patch');
    for sl = 1:2
        patchHand(sl).FaceColor = rgb(pieCols{sl});
    end
    
    legend(labels, 'Location', 'southeastoutside')
    
    %     ttl = [group(g).name ' (n = ' num2str(size(ForRevAll{g},1)) ' events) '  ttlAdd{rInd}];
    %     title(ttl)
end %group

rInd = 2;
for g = 1:2
    subplot(2,2,spMap(rInd,g))
    
    sumEv = sum(ForRevR2{g});
    if g == 1
        contData = sumEv; %control
    else
        koData = sumEv; %KO
    end %for chi squ test
    h = pie([sumEv(1) sumEv(2)]);
    
    patchHand = findobj(h, 'Type', 'Patch');
    for sl = 1:2
        patchHand(sl).FaceColor = rgb(pieCols{sl});
    end
    %     ttl = [group(g).name ' (n = ' num2str(size(ForRevR2{g},1)) ' events) '  ttlAdd{rInd}];
    %     title(ttl)
end %group

if saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    cd(curDir)
end

%% FIG 4 - SLOPES - NOT RESCALED

figtitle = 'ReplayEvents_Slopes';

if downSampEvents == 1 %down sample both rescaled and original
    minEv = min([length(slopesAll{1}) length(slopesAll{2})]);
    
    for g = 1:2
        if length(slopesAll{g}) ~= minEv
            y = datasample(rsE,1:length(slopesAll{g}),minEv,'Replace',false);
            y = sort(y);
            slopesAll{g} = slopesAll{g}(y);
            normSlopesAll{g} = normSlopesAll{g}(y);
        end %if this is the one we down sample
    end %group
    
    for sInd = 1:2
        minEv = min([length(slopesbyEvType{1,sInd}) length(slopesbyEvType{2,sInd})]);
        
        for g = 1:2
            if length(slopesbyEvType{g,sInd}) ~= minEv
                y = datasample(rsE,1:length(slopesbyEvType{g,sInd}),minEv,'Replace',false);
                y = sort(y);
                slopesbyEvType{g,sInd} = slopesbyEvType{g,sInd}(y);
            end %if this is the one we down sample
        end %group
    end %sleepInd
end %down sample events

if downSampCell == 1
    figtitle = [figtitle '_downSamp_' num2str(newCellNum) 'cell'];
end
if downSampEvents == 1
    figtitle = [figtitle '_downSampled_events'];
end

figure('Name', figtitle, 'Position', [521 191 839 746])
fmr1CircTrack_x_makeReplayEventPlots(slopesAll, slopesbyEvType)

for sp = 1:4
    subplot(2,2,sp)
    ylabel(['Slope (|deg/s|)'])
end

if saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    cd(curDir)
end

%% FIG 5 - SLOPES - RESCALED

figtitle = 'ReplayEvents_SlopesRescaled';

if downSampEvents == 1 %down sample both rescaled and original
    minEv = min([length(normSlopesAll{1}) length(normSlopesAll{2})]);
    figtitle = [figtitle '_downSampled_events'];
    
    for g = 1:2
        if length(normSlopesAll{g}) ~= minEv
            y = datasample(rsE,1:length(normSlopesAll{g}),minEv,'Replace',false);
            y = sort(y);
            normSlopesAll{g} = normSlopesAll{g}(y);
        end %if this is the one we down sample
    end %group
    
    for sInd = 1:2
        minEv = min([length(normSlopesbyEvType{1,sInd}) length(normSlopesbyEvType{2,sInd})]);
        for g = 1:2
            if length(normSlopesbyEvType{g,sInd}) ~= minEv
                y = datasample(rsE,1:length(normSlopesbyEvType{g,sInd}),minEv,'Replace',false);
                y = sort(y);
                normSlopesbyEvType{g,sInd} = normSlopesbyEvType{g,sInd}(y);
            end %if this is the one we down sample
        end %group
    end %slope Ind
end %down sample events

figure('Name', figtitle, 'Position', [521 191 839 746])
fmr1CircTrack_x_makeReplayEventPlots(normSlopesAll, normSlopesbyEvType)

for sp = 1:4
    subplot(2,2,sp)
    ylabel(['Normanlized slope (|deg/event|)'])
end

if saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    cd(curDir)
end

%% FIG 6 - DURATION

figtitle = 'ReplayEvents_duration';

if downSampEvents == 1
    minEv = min([length(eventDurs{1}) length(eventDurs{2})]);
    
    for g = 1:2
        if length(eventDurs{g}) ~= minEv
            y = datasample(rsE,1:length(eventDurs{g}),minEv,'Replace',false);
            y = sort(y);
            eventDurs{g} = eventDurs{g}(y);
        end %if we need to down samp
    end %group
    
    for sInd = 1:2
        minEv = min([length(evDursbyEvType{1,sInd}) length(evDursbyEvType{2,sInd})]);
        for g = 1:2
            if length(evDursbyEvType{g,sInd}) ~= minEv
                y = datasample(rsE,1:length(evDursbyEvType{g,sInd}),minEv,'Replace',false);
                y = sort(y);
                evDursbyEvType{g,sInd} = evDursbyEvType{g,sInd}(y);
            end %if this is the one we down sample
        end %group
    end %slope Ind
end %if down samp

if downSampCell == 1
    figtitle = [figtitle '_downSamp_' num2str(newCellNum) 'cell'];
end
if downSampEvents == 1
    figtitle = [figtitle '_downSampled_events'];
end

figure('Name', figtitle, 'Position', [521 191 839 746])
fmr1CircTrack_x_makeReplayEventPlots(eventDurs, evDursbyEvType)
for sp = 1:4
    subplot(2,2,sp)
    ylabel(['Duration (s)'])
end

if saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    cd(curDir)
end

%% PATH DISTANCE

figtitle = 'ReplayEvents_pathDistance';

if downSampEvents == 1
    minEv = min([length(pathDistAll{1}) length(pathDistAll{2})]);
    
    for g = 1:2
        if length(pathDistAll{g}) ~= minEv
            y = datasample(rsE,1:length(pathDistAll{g}),minEv,'Replace',false);
            y = sort(y);
            pathDistAll{g} = pathDistAll{g}(y);
        end %if we need to down samp
    end %group
    
    for sInd = 1:2
        minEv = min([length(pathDistbyEvType{1,sInd}) length(pathDistbyEvType{2,sInd})]);
        for g = 1:2
            if length(pathDistbyEvType{g,sInd}) ~= minEv
                y = datasample(rsE,1:length(pathDistbyEvType{g,sInd}),minEv,'Replace',false);
                y = sort(y);
                pathDistbyEvType{g,sInd} = pathDistbyEvType{g,sInd}(y);
            end %if this is the one we down sample
        end %group
    end %slope Ind
end %if down samp

if downSampCell == 1
    figtitle = [figtitle '_downSamp_' num2str(newCellNum) 'cell'];
end
if downSampEvents == 1
    figtitle = [figtitle '_downSampled_events'];
end

figure('Name', figtitle, 'Position', [521 191 839 746])
fmr1CircTrack_x_makeReplayEventPlots(pathDistAll, pathDistbyEvType)
for sp = 1:4
    subplot(2,2,sp)
    ylabel(['Path distance (deg)'])
end

if saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    cd(curDir)
end

%% JUMP DIST

figtitle = 'ReplayEvents_jumpDistance';

if downSampEvents == 1
    minEv = min([length(avgJumpAll{1}) length(avgJumpAll{2})]);
    
    for g = 1:2
        if length(avgJumpAll{g}) ~= minEv
            y = datasample(rsE,1:length(avgJumpAll{g}),minEv,'Replace',false);
            y = sort(y);
            avgJumpAll{g} = avgJumpAll{g}(y);
        end %if we need to down samp
    end %group
    
    for sInd = 1:2
        minEv = min([length(avgJumpbyEvType{1,sInd}) length(avgJumpbyEvType{2,sInd})]);
        for g = 1:2
            if length(avgJumpbyEvType{g,sInd}) ~= minEv
                y = datasample(rsE,1:length(avgJumpbyEvType{g,sInd}),minEv,'Replace',false);
                y = sort(y);
                avgJumpbyEvType{g,sInd} = avgJumpbyEvType{g,sInd}(y);
            end %if this is the one we down sample
        end %group
    end %slope Ind
end %if down samp

if downSampCell == 1
    figtitle = [figtitle '_downSamp_' num2str(newCellNum) 'cell'];
end
if downSampEvents == 1
    figtitle = [figtitle '_downSampled_events'];
end

figure('Name', figtitle, 'Position', [521 191 839 746])
fmr1CircTrack_x_makeReplayEventPlots(avgJumpAll, avgJumpbyEvType)
for sp = 1:4
    subplot(2,2,sp)
    ylabel(['Average spatial jump distance (deg)'])
end

if saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    cd(curDir)
end


%% FIG 7 - EVENT DUR x (norm) SLOPE

figtitle = 'ReplayEvents_DurationxSlope';
if downSampCell == 1
    figtitle = [figtitle '_downSamp_' num2str(newCellNum) 'cell'];
end
if downSampEvents == 1
    figtitle = [figtitle '_downSampled_events'];
    
    minEv = min([size(evDurxDist{1},1) size(evDurxDist{2},1)]);
    
    for g = 1:2
        if length(eventDurs{g}) ~= minEv
            y = datasample(rsE,1:size(evDurxDist{1},1),minEv,'Replace',false);
            y = sort(y);
            evDurxDist{g} = evDurxDist{g}(y,:);
        end %if we need to down samp
    end %group
end

figure('Name', figtitle, 'Position', [419 457 1078 420])

for g = 1:2
    subplot(1,2,g)
    
    plot(evDurxDist{g}(:,1), evDurxDist{g}(:,2), '.', 'Color', rgb(cols{g}))
    fitLine = polyfit(evDurxDist{g}(:,1), evDurxDist{g}(:,2) ,1);
    
    hold on;
    xReg = [min(evDurxDist{g}(:,1)) max(evDurxDist{g}(:,1))];
    yReg = xReg .* fitLine(1) + fitLine(2);
    plot(xReg, yReg, 'Color', rgb(cols{g}))
    
    %get r2
    
    SSR = 0;
    SSTO = 0;
    
    sampMean = mean(evDurxDist{g}(:,2));
    for i = 1:size(evDurxDist{g},1)
        x = evDurxDist{g}(i,1);
        y = evDurxDist{g}(i,2);
        
        esty = x .* fitLine(1) + fitLine(2);
        
        SSR = SSR + (esty - sampMean)^2;
        SSTO = SSTO + (y - sampMean)^2;
        
    end %i
    
    tmpRsqu = SSR / SSTO;
    
    ttl = [group(g).name ' (n = ' num2str(size(evDurxDist{g},1)) '): slope = ' num2str(fitLine(1)) ' R^2 = ' num2str(tmpRsqu)];
    title(ttl)
    
    xlabel('Event duration (s)')
    ylabel('Path distance (deg)')
    
    xlim([0 0.5])
    ylim([0 360])
    yticks(0:90:360)
    
end %group

if saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    cd(curDir)
end

%% FIG 7 - EXAMPLE PLOTS - HIGHEST
if plotHighest == 1
    if downSampCell == 0
        figtitle = 'ReplayFidelity_examples_highest';
    else
        figtitle = ['ReplayFidelity_examples_highest_downSamp_' num2str(newCellNum) 'cell'];
    end
    figure('Name', figtitle, 'Position', [248 448 1349 532])
    spMap = [1:5; 6:10];
    
    for g = 1:2
        
        tmpPpmForPlot = ppmForPlot(:,g);
        tmpR2ForPlot = r2ForPlot(:,g);
        
        [~, sortOrd] = sort(tmpR2ForPlot, 'descend');
        
        for i = 1:5
            subplot(2,5,spMap(g,find(sortOrd == i)))
            
            maxTm = size(tmpPpmForPlot{i},2) * bayesStep + bayesWin/2;
            
            imagesc(0:1, 0:360, tmpPpmForPlot{i})
            axis xy
            axis square
            colormap hot
            caxis([0 0.25])
            
            xticks([0 1])
            xticklabels({'', num2str(maxTm)})
            if g == 2
                xlabel('Time (s)')
            end
            
            yticks([0 180 360])
            if find(sortOrd == i) == 1
                ylabel({group(g).name, 'Angular position (deg)'})
            end
            
            title(['r^2 = ' num2str(round(tmpR2ForPlot(i),2))])
            
        end %which map

    end %group
            cbr = colorbar;
        set(cbr, 'Position', [.92 .60 .02 .32])
        ylabel(cbr, 'Probability')
    
    if saveOrNot == 1
        cd(saveDir)
        saveas(gcf, figtitle, 'epsc');
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        cd(curDir)
    end
end %plot highest

%% FIG 8 - EXAMPLE PLOTS - RANDOM ABOVE 0.7

if plotRandom == 1
    if downSampCell == 0
        figtitle = 'ReplayFidelity_examples_random';
    else
        figtitle = ['ReplayFidelity_examples_random_downSamp_' num2str(newCellNum) 'cell'];
    end
    figure('Name', figtitle, 'Position', [680 78 708 900])
    spMap = [1:2:9; 2:2:10]';
    
    for g = 1:2
        
        tmpPpmForPlot = ppmForPlotRand(:,g);
        tmpR2ForPlot = r2ForPlotRand(:,g);
        
        [~, sortOrd] = sort(tmpR2ForPlot, 'descend');
        
        for i = 1:5
            subplot(5,2,spMap(find(sortOrd == i),g))
            
            maxTm = size(tmpPpmForPlot{i},2) * bayesStep + bayesWin/2;
            
            imagesc(0:1, 0:360, tmpPpmForPlot{i})
            axis xy
            axis square
            colormap hot
            caxis([0 0.25])
            
            %             if max(ppmForPlot{i,g}) > 0.2
            %                 keyboard
            %             end
            
            xticks([0 1])
            xticklabels({'', num2str(maxTm)})
            if find(sortOrd == i) == 5
                xlabel('Time (s)')
            end
            
            yticks([0 180 360])
            
            if g == 1
                ylabel('Angular position (rad)')
            end
            
            if find(sortOrd == i) == 1
                title({group(g).name, ['r^2 = ' num2str(round(tmpR2ForPlot(i),2))]})
            else
                title(['r^2 = ' num2str(round(tmpR2ForPlot(i),2))])
            end
        end %subplot for group
        
        cbr = colorbar;
        set(cbr, 'Position', [.8314 .8 .03 .1])
        ylabel(cbr, 'Probability')
        
    end %group
    
    if saveOrNot == 1
        cd(saveDir)
        saveas(gcf, figtitle, 'epsc');
        saveas(gcf, figtitle, 'fig');
        saveas(gcf, figtitle, 'png');
        cd(curDir)
    end
end %plot random

%% STATS?
% preparing output for SPSS

if prepForStats == 1
    
    %     slopesForStats = zeros(length(slopesAll{1})+length(slopesAll{2}),2);
    statsSlopes = [];
    statsRsSlopes = [];
    statRsSlopesbySlope = [];
    statSlopesBySlope = [];
    statPathDistBySlope = [];
    statAvgJumpbySlope = [];
    
    statEventDurs = [];
    statEvDurbySlope = [];
    
    statForRevAll = [];
    
    statSlopeDirByGroup = [];
    
    for g = 1:2
        for i = 1:length(slopesAll{g})
            statsSlopes = [statsSlopes; g slopesAll{g}(i)];
            statsRsSlopes = [statsRsSlopes; g normSlopesAll{g}(i)];
            
            statSlopeDirByGroup = [statSlopeDirByGroup; g find(ForRevR2{g}(i,:))];
            
        end %i - slopes
        
        for i = 1:length(eventDurs{g})
            statEventDurs = [statEventDurs; g eventDurs{g}(i)];
        end %i - event Durs
        
        for sInd = 1:2
            for i = 1:length(slopesbyEvType{g,sInd})
                statSlopesBySlope = [statSlopesBySlope; g sInd slopesbyEvType{g,sInd}(i)];
                statRsSlopesbySlope = [statRsSlopesbySlope; g sInd normSlopesbyEvType{g,sInd}(i)];
                statPathDistBySlope = [statPathDistBySlope; g sInd pathDistbyEvType{g,sInd}(i)];
                statAvgJumpbySlope = [statAvgJumpbySlope; g sInd avgJumpbyEvType{g,sInd}(i)];
                
                statEvDurbySlope = [statEvDurbySlope; g sInd evDursbyEvType{g,sInd}(i)];
            end %i - slopes by slope
        end %sInd
        
    end %group
    keyboard
end %prep for stats


end %function