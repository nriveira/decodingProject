function fmr1CircTrack_x_plotSWRExample(group)
% function fmr1CircTrack_x_plotSWRExample(group)
% 
% PURPOSE:
%   This function produces figures showing sample SWR/replay events, to be
%   used as examples for illustrative purposes.
% 
% INPUT:
%   group = data struct
% 
% OUTPUT:
%   Figures. Option to save.
% 
% OPTIONS:
%   figsPerDay = how many events to be plotted for each day.
%   saveOrNot = whether to save figs (1) or not
%   r2Thresh = min threshold for r2 value for replay event. Used so that
%       examples plotted can be clear and interpretable.
%   See function for other options.
% 
% MMD
% 10/2021
% Colgin Lab

%% OPTIONS

figsPerDay = 20; %number of figs to make per day
saveOrNot = 1; %to save figs

minFr = 1; %Hz, for cell to be included

minCell = 5; %min number of cells to be active in an event

r2Thresh = 0.7; %set it a little higher, because this is just for plotting examples

%% INITIALIZE

curDir = pwd; %for returning after saving figs
saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\REPLAY\SWR_examplePlots';

lfpCut = [150 250];

%% GET DATA
for g = 1:2
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            
            figTtlBase = [group(g).name '_' group(g).rat(r).name '_' group(g).rat(r).day(d).name ];
            figCntr = 0; %reset for this day
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
            
            uIDs(badU,:) = [];
            
            for s = 2:length(group(g).rat(r).day(d).sleep) %after experience, so just sleep 2 on
                
                if isempty(group(g).rat(r).day(d).sleep(s).unit) %if there is anything in the sleep
                    continue
                end %anything in sleep
                
                ripEvents = group(g).rat(r).day(d).sleep(s).rip; %shorten variable name
                
                for i = 1:length(ripEvents)
                    
                    if ripEvents(i).r2 < r2Thresh || isnan(ripEvents(i).r2)
                        continue
                    end %doesn't meet threshold
                    
                    startTm = ripEvents(i).adjTms(1);
                    endTm = ripEvents(i).adjTms(2);
                    
                    spkTmsByU = {};
                    spkTmIDs = [];
                    spkTmPkPos = [];
                    uCntr = 0;
                    
                    maxSpk = 0;
                    maxSpkTet = [];
                    for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                        uID = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID;
                        if ismember(uID, uIDs, 'row')
                            uCntr = uCntr + 1;
                            allSpkTms = group(g).rat(r).day(d).sleep(s).unit(u).spkTms;
                            tmpEvSpkTms = allSpkTms(allSpkTms>=startTm & allSpkTms<=endTm);
                            
                            if ~isempty(tmpEvSpkTms)
                                timePassed = tmpEvSpkTms - startTm;
                                spkTmsByU = cat(1, spkTmsByU, timePassed);
                                spkTmIDs = [spkTmIDs; uCntr]; %for getting ratemap later
                                
                                [~, tmpPkPos] = max(rateMaps(uCntr,:));
                                spkTmPkPos = [spkTmPkPos; tmpPkPos];
                                
                                if length(tmpEvSpkTms) > maxSpk
                                    maxSpk = length(tmpEvSpkTms);
                                    maxSpkTet = uID(1);
                                end %more than max so dar
                            end %there were any spikes from this unit
                        end %unit was kept
                    end %unit
                    if length(spkTmsByU) < minCell
                        continue
                    end
                    
                    figCntr = figCntr + 1;
                    if figCntr > figsPerDay
                        break
                    end %made enough figs for today
                    
                    pxn = ripEvents(i).pxn;
                    
                    cd(group(g).rat(r).day(d).sleep(s).dir)
                    lfpStruct = read_in_lfp(['CSC' num2str(maxSpkTet) '.ncs']);
                    lfpStruct.data =  detrend(lfpStruct.data);
                    
                    startInd = find(lfpStruct.ts <= startTm, 1, 'Last');
                    endInd = find(lfpStruct.ts <= endTm, 1, 'Last');
                    evLFP = lfpStruct.data(startInd:endInd);
                    
                    [B,A]=butter(3,[lfpCut(1)/(lfpStruct.Fs*0.5) lfpCut(2)/(lfpStruct.Fs*0.5)]);
                    filtLfp = filtfilt(B, A, evLFP);
                    
                    figtitle = [figTtlBase '_' num2str(figCntr)];
                    figure('Name', figtitle, 'Position', [680 137 736 841])
                    
                    subplot(3,2,1)
                    plot(0:1/lfpStruct.Fs:endTm-startTm, evLFP(1:length(0:1/lfpStruct.Fs:endTm-startTm)), 'LineWidth', 1.5)
                    xlim([0 endTm-startTm])
                    xlabel('Time (s)')
                    ylabel('Voltage (mV)')
                    title('Raw LFP')
                    
                    subplot(3,2,2)
                    plot(0:1/lfpStruct.Fs:endTm-startTm, filtLfp(1:length(0:1/lfpStruct.Fs:endTm-startTm)), 'LineWidth', 1.5)
                    xlim([0 endTm-startTm])
                    xlabel('Time (s)')
                    ylabel('Voltage (mV)')
                    title('Filtered LFP (150-250 Hz)')
                    
                    subplot(3,2,4)
                    
                    [~, sortOrd] = sort(spkTmPkPos);
                    cols = hsv(length(spkTmsByU));
                    
                    hold on;
                    for u = 1:length(spkTmsByU)
                        yVal = find(u == sortOrd);
                        for st = 1:length(spkTmsByU{u})
                            line([spkTmsByU{u}(st) spkTmsByU{u}(st)], [yVal-0.3 yVal+0.3], 'Color', cols(yVal,:), 'LineWidth', 1.5);
                        end %spkTm
                    end %unit
                    
                    xlim([0 endTm-startTm])
                    xlabel('Time (s)')
                    ylim([0.5 length(spkTmsByU)+0.5])
                    if length(spkTmsByU) < 15
                        yticks(1:1:length(spkTmsByU))
                    else
                        yticks(1:2:length(spkTmsByU))
                    end %how many tick marks cna we fit
                    ylabel('Cell number')
                    
                    subplot(3,2,5)
                    hold on;
                    for u = 1:length(spkTmsByU)
                        tmpMap = rateMaps(spkTmIDs(u),:);
                        yVal = find(u == sortOrd);
                        
                        plot(tmpMap, 1:4:360, 'Color', cols(yVal,:), 'LineWidth', 1.5)
                        set(gca,'XDir','reverse');
                    end %unit
                    
                    xlabel('Firing rate (Hz)')
                    yticks('')
                    ylim([0 360])
                    
                    subplot(3,2,6)
                    imagesc(0:0.01:0.01*size(pxn,2), 1:4:360, pxn)
                    axis xy
                    colormap(hot)
                    ylabel('Position (?)')
                    xticks('')
                    xlabel(['r^2 = ' num2str(round(ripEvents(i).r2,2))])
                    %                     try
                    %                     xlabel({['r^2 = ' num2str(round(ripEvents(i).r2,2))], ['slope = ' num2str(round(ripEvents(i).slope,2)) ' ?/s']})
                    %                     catch; keyboard; end
                    cbr = colorbar;
                    set(cbr, 'Position', [.92 0.12 .02 .2])
                    ylabel(cbr, 'Probability')
                    
                    if saveOrNot == 1
                        cd(saveDir)
                        saveas(gcf, figtitle, 'epsc');
                        saveas(gcf, figtitle, 'fig');
                        saveas(gcf, figtitle, 'png');
                    end %save or not
                end %i
            end %sleep
            close all
        end %day
    end %r
end %group

cd(curDir)



end %function