function fmr1CircTrack_x_plotSingleCell_replay(group)
% function fmr1CircTrack_x_plotSingleCell_replay(group)
%
% PURPOSE:
%   To plot single cell analyses for replay events from WT and KO rats.
%
% INPUTS:
%   group: data struct
%
% OUTPUTS:
%   Figures.
%
% OPTIONS:
%   saveOrNot: save figs (1) or don't (0)
%   See function for other options.
%
% MMD
% 7/2021
% Colgin Lab

%% OPTIONS

saveOrNot = 0;
saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\REPLAY\replayEvents_singleCell';
curDir = pwd;

prepForStats = 0;

downSampEvents = 0; %so number of events equal between groups

% Binned firing rate options
binSz = 0.001; %seconds, for binning firing rate
gWinStd = 10/1000; %5 ms in Ernie's paper
gWinDur = gWinStd * 6;
nboot = 10000; %for calculating confidence intervals, from Hwaun & Colgin 2019
minFr = 1; %min peak FR in Hz to be considered active during waking behaviors (as in Hwaun & Colgin 2019)


%% INITIALIZE

spatBinSz = 4;  %for ratemaps

% Binned firing rate stuff
totTime = 0.8; %seconds, to plot binned firing rate (with event in middle)
preTime = [0.4 0.1]; %time frame before ripple to get normalizing firing rate; 0.4-0.1s in Hwaun & Colgin 2017
numBins = totTime/binSz +1;
onsetBin = ceil(numBins/2);
preBins = [totTime/2 - preTime(1)+1:onsetBin - preTime(2)/binSz];

binFr = cell(2,1); %by group
binFrPreNorm = cell(2,1); %by group - firing rate normalized by firing avg firing rate 0.4-0.1 s before replay event

gWinDur = gWinDur/binSz; %convert based on bin size
gWinStd = gWinStd/binSz; %convert based on bin size
gKrnl = gausskernel(gWinDur, gWinStd);

% FR and SPK/CELL stuff
frInEv = cell(2,1);
spkPerCell = cell(2,1);

%frac participation stuff
fracPart = cell(2,1);

% pk pos x in ev firing rate stuff
frPkPos = cell(2,360/spatBinSz);
spkPerCellPkPos = cell(2,360/spatBinSz);

% Misc
cols = {'Blue', 'Red'};
evCntr = zeros(1,2);

%% PREP FOR DOWN SAMPLE

if downSampEvents == 1
    rsE = RandStream('mt19937ar');
    
    gNum = zeros(1,2);
    
    for g = 1:2
        for r = 1:length(group(g).rat)
            for d = 1:length(group(g).rat(r).day)
                for s = 2:5
                    if ~isempty(group(g).rat(r).day(d).sleep(s).unit)
                        gNum(g) = gNum(g) + length(group(g).rat(r).day(d).sleep(s).popEvents);
                    end %if there are events
                end %sleep
            end %day
        end %rat
    end %group
    
    [minEv, minG] = min(gNum);
    keepEvs = datasample(rsE, 1:max(gNum), minEv, 'Replace', false);
    keepEvs = sort(keepEvs);
    
    allEvCntr = 0; %initialize
    
end %if down sampling events

%% GET DATA

for g = 1:2
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            
            rewLoc = group(g).rat(r).day(d).rewLocs(1);
            [~,rewInd] = min(abs(circ_dist(deg2rad(group(g).rat(r).day(d).binCtrs), deg2rad(rewLoc))-0));
            shiftVal = 1 - rewInd;
            
            numEvents = 0;
            if downSampEvents == 1 && g ~= minG %find out which events from this day to keep
                
                for s = 2:5 %getting getting total number of events for day
                    if ~isempty(group(g).rat(r).day(d).sleep(s).unit)
                        numEvents = numEvents + size(group(g).rat(r).day(d).sleep(s).popEvents,1);
                    end %if there is anything in popEvents
                end %sleep
                
                dayEvs = 1:numEvents; %initialize
                keepInds = find(keepEvs > allEvCntr & keepEvs <= allEvCntr+numEvents); %which from keepEvs are in this day (no = in first since allEvCntr starts at 0)
                dayEvInds = keepEvs(keepInds) - allEvCntr;
                dayEvs = dayEvs(dayEvInds); %the events from this day that get to be included
                
                allEvCntr = allEvCntr + numEvents;
            end %down sample
            
            countEv = 1; %so we count the events on the first pass
            for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                if max(group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap) < minFr
                    continue
                end %cell doesn't reach min fr
                [~, uPkPos] = max(group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap);
                shiftPkPos = uPkPos + shiftVal;
                if shiftPkPos < 0
                    shiftPkPos = shiftPkPos + 360/spatBinSz;
                end %if cross 0
                
                uBinFr = []; %initialize
                
                uSWRFr = []; %initialize
                uSpkPerRip = [];
                dayEvCntr = 0;
                
                uPart = 0; %for counting # events participated in
                totDen = 0; %for counting total events for the denominator of fractional participation for this unit
                for s = 2:5 %start from 2, after experience
                    if ~isempty(group(g).rat(r).day(d).sleep(s).unit) %if there is anything in the sleep
                        popEvents = group(g).rat(r).day(d).sleep(s).popEv;
                        spkTms = group(g).rat(r).day(d).sleep(s).unit(u).spkTms;
                        
                        for i = 1:size(popEvents,1)
                            dayEvCntr = dayEvCntr + 1;
                            if downSampEvents == 0 || g == minG || g ~= minG && ismember(dayEvCntr, dayEvs)
                                totDen = totDen + 1;
                                
                                startInt = popEvents(i).tms(1) - totTime/2; %start time of interest
                                endInt = popEvents(i).tms(1) + totTime/2; %end time of interest
                                
                                dur = popEvents(i).tms(2) - popEvents(i).tms(1);
                                spksSWR = find(spkTms >= popEvents(i).tms(1) & spkTms <= popEvents(i).tms(2));
                                tmpFr = length(spksSWR)/dur;
                                
                                if spksSWR ~= 0 %if cell was active in this event
                                    uSWRFr = [uSWRFr tmpFr];
                                    uSpkPerRip = [uSpkPerRip length(spksSWR)];
                                    
                                    uPart = uPart + 1;
                                end
                                
                                if startInt > 0 && endInt < group(g).rat(r).day(d).sleep(s).coords(end,1)
                                    edges = startInt:binSz:endInt; %get bin edges
                                    
                                    tmpBinCnts = histcounts(spkTms, edges);
                                    tmpBinFr = tmpBinCnts * (1/binSz); %convert to Hz
                                    uBinFr = [uBinFr; tmpBinFr]; %store for this unit
                                    
                                end %if the window fits in this begin
                                
                                if countEv == 1 %if this is the first unit (so events don't get counted for every iteration of these events)
                                    evCntr(g) = evCntr(g) + 1;
                                    %                                     if evCntr(g) > minEv
                                    %                                         keyboard
                                    %                                     end %check
                                end %whether to count event
                                
                            end %which to include (depending on down sampling events or not)
                        end %pop events
                        
                    end %if there are popevents
                end %sleep
                
                if totDen ~= 0
                    fracPart{g} = [fracPart{g} uPart/totDen];
                end
                
                if ~isempty(uBinFr)
                    binFr{g} = [binFr{g}; mean(uBinFr,1)];
                    preEvFr = mean(mean(uBinFr(:,preBins)));
                    binFrPreNorm{g} = [binFrPreNorm{g}; mean(uBinFr,1)/preEvFr*100]; %convert to %
                end %bin firing rate wasn't empty
                
                if ~isnan(mean(uSWRFr)) %as long as there were some events that had spikes for this unit, to get a mean
                    frInEv{g} = [frInEv{g} mean(uSWRFr)];
                    spkPerCell{g} = [spkPerCell{g} mean(uSpkPerRip)];
                    try
                    frPkPos{g,shiftPkPos} = [frPkPos{g,shiftPkPos} mean(uSWRFr)];
                    catch
                        keyboard
                    end
                    spkPerCellPkPos{g,shiftPkPos} = [spkPerCellPkPos{g,shiftPkPos} mean(uSpkPerRip)];
                end %mean isn't nan
                countEv = 0;
            end %unit
        end %day
    end %rat
end %group

%% FIG 1 - BINNED FR AROUND RIP
keyboard
figtitle = 'BinnedFr_inReplayEvent';

if downSampEvents == 1
    figtitle = [figtitle '_downSampled_events'];
end %down samp
figure('Position', [333 501 1142 388], 'Name', figtitle)

yLabs = {'Firing rate (Hz)', 'Normalized firing rate (%)'};

for sp = 1:2
    subplot(1,2,sp)
    lh = zeros(1,2);
    leg = cell(1,2);
    
    for g = 1:2
        if sp == 1
            pullData = binFr{g};
        else
            pullData = binFrPreNorm{g};
        end %which subplot
        cellNum = size(pullData,1);
        
        mubin = mean(pullData);
        smbin = conv(mubin, gKrnl, 'same'); %smooth it
        smbinCut = smbin(onsetBin-(totTime/2/binSz):onsetBin+(totTime/2/binSz)-1);
        
        CI = bootci(nboot, {@mean, pullData}, 'type', 'per');
        for row = 1:2
            smCI(row,:) = conv(CI(row,:), gKrnl, 'same');
        end
        smCIcut = smCI(:,onsetBin-(totTime/2/binSz):onsetBin+(totTime/2/binSz)-1);
        
        hold on;
        lh(g) =  plot_filled_ci(1:(totTime/binSz), smbinCut, smCI, rgb(cols{g}));
        
        leg{g} = [group(g).name ' n = ' num2str(cellNum) ' cells\newline(' num2str(evCntr(g)) ' events)'];
    end %group
    
    ylabel(yLabs{sp})
    xlabel('Time from replay event onset (s)')
    xlim([0 numBins])
    % ylim([0 5.5])
    xticks([0 ((totTime/binSz)+1)/2 (totTime/binSz)+1])
    xticklabels(-totTime/2:totTime/2:totTime/2)
    
    legend(lh, leg, 'Location', 'northwest')
end %subplot

if saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    cd(curDir)
end

%% FIG 2 - FR IN RIP

figtitle = 'FiringRate_inReplayEvent';
if downSampEvents == 1
    figtitle = [figtitle '_downSampled_events'];
end %down samp
figure('Name', figtitle, 'Position', [574 538 790 368]);

findMax = max([max(frInEv{1}) max(frInEv{2})]);
findMin = min([min(frInEv{1}) min(frInEv{2})]);

legText = {group(1).name, group(2).name};

subplot(1,2,1)
for g = 1:2
    minFr = floor(findMin);
    maxFr = ceil(findMax) + 3;
    
    sigma = 2;
    
    [Wdis_pre,axisValue1] = WeightedProportion(frInEv{g}', minFr, maxFr, sigma);
    %     if min(frInEv{g}) ~= findMin
    plot(axisValue1, Wdis_pre, 'Color', rgb(cols{g}), 'LineWidth', 1)
    %     else
    %         plot([floor(findMin) axisValue1], [0 Wdis_pre], 'Color', rgb(cols{g}), 'LineWidth', 1)
    %     end
    
    
    legText{g} = [legText{g} ' n = ' num2str(length(frInEv{g,1})) ' cells\newline(' num2str(evCntr(g)) ' events)'];
    
    hold on;
    
end %group
ylim([0 0.06])
xticks(0:10:30)

ax = gca;
yMax = ax.YLim(2);

medFr = cellfun(@median, frInEv);
for g = 1:2
    line([medFr(g) medFr(g)], [0 yMax], 'Color', rgb(cols{g}), 'LineStyle', '--')
end

legend(legText, 'Location', 'southoutside')
xlim([0 ceil(findMax)+3])
xlabel('Firing rate (Hz)')
ylabel('Proportion of neurons')
title('Active CA1 neurons')

subplot(1,2,2)

jitter = 0.2;
h = dotplot(1:2, frInEv, jitter, [rgb(cols{1}); rgb(cols{2})], [rgb('Black'); rgb('Black')]);

xticklabels({group(1).name, group(2).name})
ylabel('Firing rate (Hz)')

% meanFr = cellfun(@mean, frInEv);
for g = 1:2
    line([g-0.15 g+0.15], [medFr(g) medFr(g)], 'Color', 'Black', 'LineWidth', 3)
end

P = ranksum(frInEv{1}, frInEv{2});

if P < 0.05
    xlabel({['Mann-Whitney U Test: p = ' num2str(P)]})
end

if saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    cd(curDir)
end

%% FIG 3 - SPK PER SWR, PER CELL

figtitle = 'SpkPerSWR_perCell_inReplayEvent';
if downSampEvents == 1
    figtitle = [figtitle '_downSampled_events'];
end %down samp
figure('Name', figtitle, 'Position', [574 538 790 368]);

findMax = max([max(spkPerCell{1}) max(spkPerCell{2})]);
tmp = min([min(spkPerCell{1}) min(spkPerCell{2})]);

legText = {group(1).name, group(2).name};

subplot(1,2,1)
for g = 1:2
    minSpk = 1;
    maxSpk = ceil(findMax);
    
    sigma = 0.3;
    
    [Wdis_pre,axisValue1] = WeightedProportion(spkPerCell{g}', minSpk, maxSpk, sigma);
    plot(axisValue1, Wdis_pre, 'Color', rgb(cols{g}), 'LineWidth', 1) %because before 1, prop = 0
    %     plot([1 axisValue1], [0 Wdis_pre], 'Color', rgb(cols{g}), 'LineWidth', 1) %because before 1, prop = 0
    
    legText{g} = [legText{g} ' n = ' num2str(length(spkPerCell{g,1})) ' cells\newline(' num2str(evCntr(g)) ' events)'];
    hold on;
end %group

ylim([0 0.06])
xticks(0:1:4)

ax = gca;
yMax = ax.YLim(2);

medSpk = cellfun(@median, spkPerCell);
for g = 1:2
    line([medSpk(g) medSpk(g)], [0 yMax], 'Color', rgb(cols{g}), 'LineStyle', '--')
end

legend(legText, 'Location', 'southoutside')
xlim([0 ceil(findMax)])
% ylim([0 0.05])
xlabel('Spike/replay event, per cell')
ylabel('Proportion of neurons')
title('Active CA1 neurons')

subplot(1,2,2)

jitter = 0.2;
h = dotplot(1:2, spkPerCell, jitter, [rgb(cols{1}); rgb(cols{2})], [rgb('Black'); rgb('Black')]);

xticklabels({group(1).name, group(2).name})
ylabel('Spike/replay event, per cell')

% medSpk = cellfun(@median, spkPerCell);
for g = 1:2
    line([g-0.15 g+0.15], [medSpk(g) medSpk(g)], 'Color', 'Black', 'LineWidth', 3)
end

P = ranksum(spkPerCell{1}, spkPerCell{2});

if P < 0.05
    xlabel({['Mann-Whitney U Test: p = ' num2str(P)]})
end

if saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    cd(curDir)
end

%% FIG 4 - FRACTIONAL PARTICIPATION

figtitle = 'FractionalParticipation_inReplayEvent';
if downSampEvents == 1
    figtitle = [figtitle '_downSampled_events'];
end %down samp
figure('Name', figtitle, 'Position', [574 538 790 368]);

subplot(1,2,1)

jitter = 0.2;
h = dotplot(1:2, fracPart, jitter, [rgb(cols{1}); rgb(cols{2})], [rgb('Black'); rgb('Black')]);

avgPart = cellfun(@mean, fracPart);
for g = 1:2
    line([g-0.15 g+0.15], [avgPart(g) avgPart(g)], 'Color', 'Black', 'LineWidth', 3)
end

xticklabels({group(1).name, group(2).name})
ylabel('Fractional participation in replay event')

subplot(1,2,2)

bgraph = bar(avgPart, 'FaceColor', 'Flat');

semPart = cellfun(@semfunct, fracPart);
hold on;
er = errorbar(avgPart, semPart);

er.Color = rgb('Black');
er.LineStyle = 'None';
for g = 1:2
    bgraph.CData(g,:) = rgb(cols{g});
end %group

xticklabels({group(1).name, group(2).name})
ylabel('Fractional participation in replay event')

if saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    cd(curDir)
end

%% FIG 5 - FR AND SPC BY PEAK POS

figtitle = 'FiringRate_spikePerCell_byPeakPos_inReplayEvent';
if downSampEvents == 1
    figtitle = [figtitle '_downSampled_events'];
end %down samp
figure('Name', figtitle, 'Position', [574 538 790 368]);

yLabs = {'Firing rate', 'Spike/replay event, per cell'};
for sp = 1:2
    subplot(1,2,sp)
    
    if sp == 1
        pullData = frPkPos;
    else
        pullData = spkPerCellPkPos;
    end %which data we need for this subplot
    
    for g = 1:2
        meanData = cellfun(@mean, pullData(g,:));

        plot(1:spatBinSz:360, meanData, 'Color', rgb(cols{g}))
        hold on;
    end %group
    
    xlim([0 360])
    xlabel('Position (?)')
    ylabel(yLabs{sp})
    if sp == 2
        ylim([0 6])
    end
    
end %subplot


%% STATS? - set up to plug into SPSS

if prepForStats == 1
    
    STATbinFr = []; %time after replay event onset, as in Hwaun & Colgin 2019
    STATfrInEv = [];
    STATspkPerCell = [];
    STATfracPart = [];
    
    semFr = cellfun(@semfunct, frInEv);
    semSpk = cellfun(@semfunct, spkPerCell);
    
    for g = 1:2
        fprintf('%s median Fr in event = %d, sem = %d\n', group(g).name, medFr(g), semFr(g))
        fprintf('%s median spk per cell = %d, sem = %d\n', group(g).name, medSpk(g), semSpk(g))
        
        for i = 1:size(binFr{g},1) %same for all - same cells
            STATbinFr = [STATbinFr; g binFr{g}(i,onsetBin:end)];
            
            STATfrInEv = [STATfrInEv; g frInEv{g}(i)];
            STATspkPerCell = [STATspkPerCell; g spkPerCell{g}(i)];
        end %i - cells in binFr for group
        for i = 1:length(fracPart{g})
            STATfracPart = [STATfracPart; g fracPart{g}(i)];
        end
    end %group
    
    keyboard
end %stat prep

end %function