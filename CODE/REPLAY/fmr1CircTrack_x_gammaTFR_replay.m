function fmr1CircTrack_x_gammaTFR_replay(group)
% function fmr1CircTrack_x_gammaTFR_replay(group)
%
% PURPOSE:
%   Plot TFR before and after the onset of replay event for WT and KO rats.
%   See function for citation/clear description of method. Also, this code
%   can take a long time.
%
% INPUT:
%   group struct
%
% OUTPUT:
%   Figures. Can be saved.
%
% OPTIONS:
%   See function for options.
%
% MMD
% 7/2021
% Colgin Lab

% NOTE: This function uses the Mably et al. 2017 method (modified). TFR of
% power for multiple replay events recorded from the same tetrode and DAY
% (not session) within the same animal were then averaged. Power esimated
% from (ALL SLEEPS) of each recording day for each tetrode with place cells
% were then averaged and z-scored. z-scored power estimates were then
% averaged across recording days for each rat. The grand average TFR of
% power produced were created by averaging across all rats within each
% genotype.

%% OPTIONS

saveOrNot = 0;
saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\REPLAY\TFRreplayEvent';

downSampEvents = 0; %to be equal between genotypes - events that are discarded will be random, not defined by animal day etc.

preEvTm = 0.1; %time pre-event
postEvTm = 0.25;

width = 7; %as in Mably et al. 2017
fRange = 25:1:95; %slow gamma

%% INITALIZE

allzTFR = cell(2,1); %by group, zscored

groupEvents = zeros(2,1);

Fs = 2000;

voltConv = 1000000; %for converting microvolts to volts

order = 5; %prewhiten

curDir = pwd;

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
    fprintf('%s\n', group(g).name)
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d\n', r, length(group(g).rat))
        
        ratzTFR = [];
        
        for d = 1:length(group(g).rat(r).day)
            
            numEvents = 0;
            for s = 2:5 %getting getting total number of events for day
                if ~isempty(group(g).rat(r).day(d).sleep(s).unit)
                    
                    numEvents = numEvents + size(group(g).rat(r).day(d).sleep(s).popEvents,1);
                end %if there is anything in popEvents
            end %sleep
            
            if downSampEvents == 1 && g ~= minG %find out which events from this day to keep
                
                dayEvs = 1:numEvents; %initialize
                keepInds = find(keepEvs > allEvCntr & keepEvs <= allEvCntr+numEvents); %which from keepEvs are in this day (no = in first since allEvCntr starts at 0)
                dayEvInds = keepEvs(keepInds) - allEvCntr;
                dayEvs = dayEvs(dayEvInds); %the events from this day that get to be included
                
                allEvCntr = allEvCntr + numEvents;
                
            end %down sample
            
            if downSampEvents == 0 || g == minG
                fprintf('\t\tDay %d/%d - %d events\n', d, length(group(g).rat(r).day), numEvents)
            else
                fprintf('\t\tDay %d/%d - %d events (downsampled)\n', d, length(group(g).rat(r).day), length(dayEvs))
            end
            
            if numEvents > 0
                groupEvents(g) = groupEvents(g) + numEvents;
                
                tetNums = group(g).rat(r).day(d).tetNums; %all tets with cells
                dayTFR = zeros(length(fRange), (postEvTm+preEvTm)*Fs+1, length(tetNums));
                
                for tt = 1:length(tetNums)
                    fprintf('\t\t\tTetrode %d/%d\n', tt, length(tetNums))
                    
                    cutBinnedLfp = zeros((postEvTm+preEvTm)*Fs+1, numEvents);
                    iCntr = 0; %for filling in events
                    
                    for s = 2:5
                        if ~isempty(group(g).rat(r).day(d).sleep(s).unit)
                            cd(group(g).rat(r).day(d).sleep(s).dir)
                            
                            lfpStruct = read_in_lfp(['CSC' num2str(tetNums(tt)) '.ncs']);
                            voltData = lfpStruct.data * voltConv; %convert to volts - idk Alex does this
                            
                            popEvents = group(g).rat(r).day(d).sleep(s).popEvents; %shorten variable name
                            
                            for i = 1:length(popEvents)
                                iCntr = iCntr + 1; %events in day
                                
                                if downSampEvents == 0 || g == minG || ismember(iCntr, dayEvs)
                                    startTm = popEvents(i,1) - preEvTm;
                                    
                                    startInd = find(lfpStruct.ts <= startTm, 1, 'Last');
                                    endInd = startInd + (preEvTm+postEvTm)*Fs;
                                    
                                    tmpLfp = voltData(startInd:endInd);
                                    
                                    whiteTmpLfp = prewhitening(tmpLfp, order);
                                    cutBinnedLfp(:,iCntr) = whiteTmpLfp;
                                    
                                    if tt == 1
                                        evCntr(g) = evCntr(g) + 1;
                                    end %so events don't get exponentially counted
                                end %if this isn't being excluded due to random down sampling
                                
                            end %popevents
                            
                        end %if anything in sleep
                    end %sleep
                    fprintf('\t\t\t\tGetting tetrode TFR...\n')
                    
                    if downSampEvents == 1
                        zeroInds = find(sum(cutBinnedLfp) == 0);
                        cutBinnedLfp(:,zeroInds) = [];
                        
                    end %down samp
                    
                    [inTetTFR,timeVec] = traces2TFR(cutBinnedLfp, fRange, Fs, width);
                    
                    dayTFR(:,:,tt) = inTetTFR;
                end %tetrode
                
                fprintf('\t\t\tGetting day TFR...\n')
                dayMeanTFR = mean(dayTFR,3);
                meanTFR = mean(dayMeanTFR(:));
                stdTFR = std(dayMeanTFR(:));
                
                dayzTFR = (dayMeanTFR - meanTFR) ./ stdTFR; %get zscore
                
                ratzTFR = cat(3, ratzTFR, dayzTFR);
            end %if there are any events for this day
        end %day
        
        fprintf('\t\tGetting rat TFR...\n')
        allzTFR{g} = cat(3, allzTFR{g}, mean(ratzTFR,3));
        
    end %rat
    
end %group

keyboard
%% FIGURE

if downSampEvents == 0
    figtitle = 'TFR_gamma_replayEvents';
else
    figtitle = 'TFR_gamma_replayEvents_downSampled';
end %down samp? for fig name

figure('Name', figtitle, 'Position', [386 423 1177 420])

cAx = zeros(1,2);

for g = 1:2
    subplot(1,2,g)
    
    pullTFR = allzTFR{g};
    % pullTFR = tmpAllTFR{g};
    pullTFR = mean(pullTFR,3);
    
    timeAx = [-preEvTm:0.05:postEvTm];
    imagesc(timeAx,fRange,pullTFR)
    axis xy
    
%     xlim([-0.2 0.2])
    xticks(-preEvTm:0.1:postEvTm)
    xlabel('Time from replay event (s)')
    line([0 0], [fRange(1) fRange(end)], 'LineStyle', '--', 'Color', 'White')
    
    yticks(fRange(1):5:fRange(end))
    ylabel('Frequency (Hz)')
    
    title([group(g).name ' (n = ' num2str(evCntr(g)) ' events)'])
    
    if max(pullTFR(:)) > cAx(2)
        cAx(2) = max(pullTFR(:));
    end
    if min(pullTFR(:)) < cAx(1)
        cAx(1) = min(pullTFR(:));
    end
    
end %group

for g = 1:2
    subplot(1,2,g)
    caxis(cAx)
end

cbr = colorbar;
cbr.Position = [.93 0.1095 0.0227 .8167];
ylabel(cbr, 'Power (z-score)')

if saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    cd(curDir)
end

end %function