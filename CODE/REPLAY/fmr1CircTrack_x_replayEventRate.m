function fmr1CircTrack_x_replayEventRate(group)
% function fmr1CircTrack_x_replayEventRate(group)
%
% PURPOSE:
%   Plot replay event rate during the time in the rest pot based on various
%   normalizing factors (total time in pot, speed < 5cm/s, etc.)
%
% INPUT:
%   group struct
%
% OUTPUT:
%   Figures.
%
% MMD
% Colgin Lab
% 8/2021

%NOTE: n IS CURRENTLY DAYS, NOT RATS. CHANGE LATER

%% INITIALIZE/OPTIONS

saveOrNot = 1;
saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\REPLAY\replayRate';

curDir = pwd;

cols = {'Blue', 'Red'};
xOffset = [-0.1 0.1];

timeInPot = cell(2,5); %rat over whole time in pot
speedLimit = cell(2,5); %rate when rat is under 5 cm/s
immPeriods = cell(2,5); %rate during "immobile periods"
thetaDeltRat = cell(2,5); %rate when theta/delta ratio < 1

runThresh = 5; %cm/s -- threshold for movement

%% GET DATA

for g = 1:2
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            for s = 1:5
                if ~isempty(group(g).rat(r).day(d).sleep(s).unit)
                    numEv = size(group(g).rat(r).day(d).sleep(s).popEvents,1);
                    
                    totTm = group(g).rat(r).day(d).sleep(s).coords(end,1) - group(g).rat(r).day(d).sleep(s).coords(1,1);
                    timeInPot{g,s} = [timeInPot{g,s} numEv/totTm];
                    
                    instRs = get_runspeed(group(g).rat(r).day(d).sleep(s).coords);
                    [mwRs, mwInds] = mw_avg(instRs(:,2),15,7,1); % moving window average the runspeed to 0.5 s windows with 0.23 s steps
                    mwRs = [instRs(1:mwInds(1)-1) mwRs instRs(mwInds(end)+1:length(instRs))]; %make the length of the moving window version match
                    velFiltBnry = zeros(1,length(mwRs));
                    velFiltBnry(mwRs<runThresh) = 1; %1 = less than 5 cm/s
                    
                    chunks = bwconncomp(velFiltBnry,4); %find where multiple consecutive frames of slow speed are
                    
                    speedTime = 0;
                    for ch = 1:chunks.NumObjects
                        curChunk = chunks.PixelIdxList{ch};
                        startInd = curChunk(1,1);
                        endInd = curChunk(end,1);
                        startTm = group(g).rat(r).day(d).sleep(s).coords(startInd,1); %convert to seconds
                        endTm = group(g).rat(r).day(d).sleep(s).coords(endInd,1);
                        
                        speedTime = speedTime + (endTm-startTm); %save it for output
                    end %chunks
                    
                    speedLimit{g,s} = [speedLimit{g,s} numEv/speedTime];
                    
                    [~, timeImm] = find_immobile_periods(group(g).rat(r).day(d).sleep(s).coords);
                    immPeriods{g,s} = [immPeriods{g,s} numEv/timeImm];
                    
                    cd(group(g).rat(r).day(d).sleep(s).dir)
                    
                    thetaTet = group(g).rat(r).day(d).thetaTet;
                    load(['CSC' num2str(thetaTet) '_LfpForRipDetect.mat']);
                    
                    deltaAmpTs = abs(hilbert(lfpStruct.deltaFiltLfp));
                    thetaAmpTs = abs(hilbert(lfpStruct.thetaFiltLfp));
                    
                    thetaDeltaRatio = thetaAmpTs ./ deltaAmpTs;
                    
                    thetaDeltaRatio(thetaDeltaRatio == Inf) =  max(thetaDeltaRatio(thetaDeltaRatio~=Inf)); %get rid of infinites
                    thetaDeltaRatio(thetaDeltaRatio == 0) = 0.0001; %get rid of zeros
                    
                    avg = nanmean(thetaDeltaRatio);
                    stdev = nanstd(thetaDeltaRatio);
                    
                    thetaDeltaRatio = (thetaDeltaRatio - avg) ./ stdev;
                    
                    timeTDRat = length(find(thetaDeltaRatio < 0 )) / lfpStruct.Fs;
                    
                    thetaDeltRat{g,s} =  [thetaDeltRat{g,s} numEv/timeTDRat];
                end %if not empty
            end %sleep
        end %day
    end %rat
end %group

keyboard
%% FIGURES

figtitle = 'Rate_timeInPot';
figure('Name', figtitle)

meanTimeInPot = cellfun(@mean, timeInPot);
semTimeInPot = cellfun(@semfunct, timeInPot, 'UniformOutput', false);

checkErr = cellfun(@isempty, semTimeInPot);
semTimeInPot(find(checkErr)) = {0};
semTimeInPot = cell2mat(semTimeInPot);

nDays = cellfun(@length, timeInPot);

for g = 1:2
    hold on;
    xVals = (1:5) + xOffset(g);
    er = errorbar(xVals, meanTimeInPot(g,:), semTimeInPot(g,:), 'Color', rgb(cols{g}));
    
end %group

ylim([0 0.5])
ylabel('Replay event rate (Hz)')

xlim([0 6])
xticks(1:5)
xlabel('Sleep #')

legend(group(1).name, group(2).name)

title('Time in pot')

if saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

%~~~~~~~~~~~~~~~~

figtitle = 'Rate_speedUnder5';
figure('Name', figtitle)

meanRate = cellfun(@mean, speedLimit);
semRate = cellfun(@semfunct, speedLimit, 'UniformOutput', false);

checkErr = cellfun(@isempty, semRate);
semRate(find(checkErr)) = {0};
semRate = cell2mat(semRate);

nDays = cellfun(@length, timeInPot);

for g = 1:2
    hold on;
    xVals = (1:5) + xOffset(g);
    er = errorbar(xVals, meanRate(g,:), semRate(g,:), 'Color', rgb(cols{g}));
end %group

ylim([0 0.5])
ylabel('Replay event rate (Hz)')

xlim([0 6])
xticks(1:5)
xlabel('Sleep #')

legend(group(1).name, group(2).name)

title('Speed <5cm/s')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

%~~~~~~~~~~~~~~~~

figtitle = 'Rate_immobilePeriods';
figure('Name', figtitle)

meanRate = cellfun(@mean, immPeriods);
semRate = cellfun(@semfunct, immPeriods, 'UniformOutput', false);

checkErr = cellfun(@isempty, semRate);
semRate(find(checkErr)) = {0};
semRate = cell2mat(semRate);

nDays = cellfun(@length, timeInPot);

for g = 1:2
    hold on;
    xVals = (1:5) + xOffset(g);
    er = errorbar(xVals, meanRate(g,:), semRate(g,:), 'Color', rgb(cols{g}));
end %group

ylim([0 0.7])
ylabel('Replay event rate (Hz)')

xlim([0 6])
xticks(1:5)
xlabel('Sleep #')

legend(group(1).name, group(2).name)

title('Immobile periods')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

%~~~~~~~~~~~~~~~~

figtitle = 'Rate_thetaDeltaRatio';
figure('Name', figtitle)

meanRate = cellfun(@mean, thetaDeltRat);
semRate = cellfun(@semfunct, thetaDeltRat, 'UniformOutput', false);

checkErr = cellfun(@isempty, semRate);
semRate(find(checkErr)) = {0};
semRate = cell2mat(semRate);

nDays = cellfun(@length, timeInPot);

for g = 1:2
    hold on;
    xVals = (1:5) + xOffset(g);
    er = errorbar(xVals, meanRate(g,:), semRate(g,:), 'Color', rgb(cols{g}));
end %group

ylim([0 0.6])
ylabel('Replay event rate (Hz)')

xlim([0 6])
xticks(1:5)
xlabel('Sleep #')

legend(group(1).name, group(2).name)

title('Theta/delta ratio')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
    cd(curDir)
end %save option



end %function