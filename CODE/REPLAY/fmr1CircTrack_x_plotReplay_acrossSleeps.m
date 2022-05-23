function fmr1CircTrack_x_plotReplay_acrossSleeps(group)
% function fmr1CircTrack_x_plotReplay_acrossSleeps(group)
%
% PURPOSE:
%   Plot the replay properties across sleep sessions to visualize any
%   changes throughout the day.
% 
% INPUT:
%   group = data struct, through fxn 6
% 
% OUTPUT:
%   Figures:
%       - replay fidelity
%       - selectivity index of forward/reverse events
%       - slopes
%       - normalize slope
%       - path distance
%       - event duration
% 
% MMD
% 11/2021
% Colgin Lab

%% OPTIONS

saveOrNot = 1; %to save figs

minFr = 1; %Hz, for cell to be included

saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\REPLAY\replayProperties\acrossSleeps';

cols = {'Blue', 'Red'};
nboot = 5000; %for calculating confidence intervals

jitter = 0.2; %for dotplots

r2Thresh = 0.5;


%% INITIALIZE

% Decoding parameters
sampRate = 20000; %HzsampRate = 20000; %Hz - spike sampling rate
bayesWin = 20/1000; %20 ms time window, as in Hwaun & Colgin 2019
bayesStep = 10/1000; %10 ms time step, as in Hwaun & Colgin 2019

degBinCtrs = group(2).rat(1).day(1).binCtrs; %doesn't change across days/rats
radBinCtrs = deg2rad(degBinCtrs);

slopeNames = {'Forward', 'Reverse'};

r2All = cell(2,5);
r2ByEvType = cell(2,5,2);

slopesAll = cell(2,5); %group x sleep #
slopesByEvType = cell(2,5,2); %group x sleep # x ev type

normSlopes = cell(2,5); %rescaled so time of SWR event is 0-1
normSlopesbyEvType = cell(2,5,2);

pathDist = cell(2,5); %just dist between end and beginning of decoded path
pathDistbyEvType = cell(2,5,2);

eventDurs = cell(2,5); %for storing event durations
evDursbyEvType = cell(2,5,2); %forwards vs. reverse

eventType = cell(2,5); %forward or reverse

frInEv = cell(2,5);
spkPerCell = cell(2,5);

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
            
            for s = 1:5
                if isempty(group(g).rat(r).day(d).sleep(s).unit) %if there is anything in the sleep
                    continue %to next sleep (not break - one rat only has data for sleeps 2-4)
                end %anything in sleep
                
                popEvents = group(g).rat(r).day(d).sleep(s).popEv; %shorten variable name
                
                for i = 1:length(popEvents)
                    r2 = popEvents(i).r2;
                    r2All{g,s} = [r2All{g,s} r2];
                    
                    startTm = popEvents(i).tms(1);
                    endTm = popEvents(i).tms(2);
                    eventDurs{g,s} = [eventDurs{g,s} endTm-startTm];
                    
                    if r2 < r2Thresh || isnan(r2)
                        continue %to next event
                    end %exceeds r2 threshold
                    
                    slope = popEvents(i).slope;
                    slopesAll{g,s} = [slopesAll{g,s}; rad2deg(abs(slope))];
                    
                    normSlope = slope .* (endTm-startTm); %rescale slope
                    normSlopes{g,s} = [normSlopes{g,s} rad2deg(abs(normSlope))];
                    
                    tmpForRev = zeros(1,2);
                    if slope > 0 %forward event
                        sInd = 1; %positive slope
                        tmpForRev(1) = 1;
                        
                    else %negative event
                        sInd = 2;
                        tmpForRev(2) = 1;
                    end %which direction is slope
                    
                    pxn = popEvents(i).pxn;
                    [maxVals, maxInds] = max(pxn);
                    notNan = find(~isnan(maxVals));
                    tmpDist = circ_dist(radBinCtrs(maxInds(notNan(end))), radBinCtrs(maxInds(notNan(1)))); %distance between two points
                    
                    if sInd == 1 && tmpDist < 0  %circ_dist just takes the shortest abs val amount, + or -
                        tmpDist = 2*pi + tmpDist;
                    elseif sInd == 2 && tmpDist > 0
                        tmpDist = 2*pi - tmpDist;
                    end %checking disntance is correct
                    pathDist{g,s} = [pathDist{g,s} rad2deg(abs(tmpDist))];
                    
                    r2ByEvType{g,s,sInd} = [r2ByEvType{g,s,sInd} r2];
                    eventType{g,s} = [eventType{g,s}; tmpForRev];
                    slopesByEvType{g,s,sInd} = [slopesByEvType{g,s,sInd} rad2deg(abs(slope))];
                    evDursbyEvType{g,s,sInd} = [evDursbyEvType{g,s,sInd} endTm-startTm];
                    normSlopesbyEvType{g,s,sInd} = [normSlopesbyEvType{g,s,sInd} rad2deg(abs(normSlope))];
                    pathDistbyEvType{g,s,sInd} = [pathDistbyEvType{g,s,sInd} rad2deg(abs(tmpDist))];
                    
                end %pop events - i
            end %sleep
        end %day
    end %rat
    
end %group
keyboard
cd(saveDir)
%% FIG 1 - REPLAY FIDELITY

figtitle = 'ReplayFidelity_acrossSleeps';
figure('Position', [717 504 525 388], 'Name', figtitle)

meanR2 = cellfun(@nanmean, r2All);
semR2 = cellfun(@nansemfunct, r2All);

lh = zeros(1,2);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:5, meanR2(g,:), semR2(g,:), semR2(g,:), 'Color', rgb(cols{g}));
end %group

xlim([0 6])
xticks(1:5)
xlabel('Sleep session')

ylim([0.2 0.5])
ylabel('r^2')

title('Replay fidelity')
legend(lh, {group(1).name, group(2).name})

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 2 - FORWARD VS REVERSE PROP
keyboard
figtitle = 'ForwardReverseProportion_acrossSleeps';
figure('Position', [717 504 525 388], 'Name', figtitle)

findEmpt = cellfun(@isempty, eventType);
eventType(findEmpt) = {zeros(1,2)};

tmpFunc = @(x)((sum(x(:,1))-sum(x(:,2)))/sum(x(:)));
eventDirInd = cellfun(tmpFunc, eventType);

bgraph = bar(eventDirInd', 'FaceColor', 'Flat');
for g = 1:2
    bgraph(g).CData = rgb(cols{g});
end %group

xlabel('Sleep')

ylim([-0.5 0.5])
ylabel('<- More reverse         More forward ->')

legend({group(1).name, group(2).name}, 'Location', 'northeastoutside')
title(['Replay events (r^2 > ' num2str(r2Thresh) ')'])

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end


%% FIG 3 - SLOPES

figtitle = 'ReplaySlope_acrossSleeps';
figure('Position', [717 504 525 388], 'Name', figtitle)

meanSlope = cellfun(@mean, slopesAll);
try
    semSlope = cellfun(@semfunct, slopesAll);
catch
    semSlope = cellfun(@semfunct, slopesAll, 'UniformOutput', false);
    findEmpt = cellfun(@isempty, semSlope);
    semSlope(findEmpt) = {0};
    semSlope = cell2mat(semSlope);
end %some empty bins

lh = zeros(1,2);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:5, meanSlope(g,:), semSlope(g,:), semSlope(g,:), 'Color', rgb(cols{g}));
end %group

xlim([0 6])
xticks(1:5)
xlabel('Sleep session')

ylim([1500 3000])
ylabel('Slope (°/s)')

title(['Replay event slope (r^2 > ' num2str(r2Thresh) ')'])
legend(lh, {group(1).name, group(2).name})

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 4 - SLOPES BY EV TYPE

figtitle = 'ReplaySlope_acrossSleeps_byEvType';
figure('Position', [472 504 1011 388], 'Name', figtitle)

meanSlope = cellfun(@mean, slopesByEvType);
try
    semSlope = cellfun(@semfunct, slopesByEvType);
catch
    semSlope = cellfun(@semfunct, slopesByEvType, 'UniformOutput', false);
    findEmpt = cellfun(@isempty, semSlope);
    semSlope(findEmpt) = {0};
    semSlope = cell2mat(semSlope);
end %some empty bins

for sInd = 1:2
    subplot(1,2,sInd)
    lh = zeros(1,2);
    for g = 1:2
        hold on;
        lh(g) = errorbar(1:5, meanSlope(g,:,sInd), semSlope(g,:,sInd), semSlope(g,:,sInd), 'Color', rgb(cols{g}));
    end %group
    
    xlim([0 6])
    xticks(1:5)
    xlabel('Sleep session')
    
    ylim([1200 3500])
    ylabel('Slope (°/s)')
    
    title([slopeNames{sInd} ' event slope (r^2 > ' num2str(r2Thresh) ')'])
    legend(lh, {group(1).name, group(2).name}, 'Location', 'northeastoutside')
end %ev direction

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 5 - NORM SLOPES

figtitle = 'ReplayNormSlope_acrossSleeps';
figure('Position', [717 504 525 388], 'Name', figtitle)

meanSlope = cellfun(@mean, normSlopes);
try
    semSlope = cellfun(@semfunct, normSlopes);
catch
    semSlope = cellfun(@semfunct, normSlopes, 'UniformOutput', false);
    findEmpt = cellfun(@isempty, semSlope);
    semSlope(findEmpt) = {0};
    semSlope = cell2mat(semSlope);
end %some empty bins

lh = zeros(1,2);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:5, meanSlope(g,:), semSlope(g,:), semSlope(g,:), 'Color', rgb(cols{g}));
end %group

xlim([0 6])
xticks(1:5)
xlabel('Sleep session')

ylabel('Normalize slope (°/event)')

title(['Replay event norm slope (r^2 > ' num2str(r2Thresh) ')'])
legend(lh, {group(1).name, group(2).name})

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 6 - NORM SLOPES BY EV TYPE

figtitle = 'ReplayNormSlope_acrossSleeps_byEvType';
figure('Position', [472 504 1011 388], 'Name', figtitle)

meanSlope = cellfun(@mean, normSlopesbyEvType);
try
    semSlope = cellfun(@semfunct, normSlopesbyEvType);
catch
    semSlope = cellfun(@semfunct, normSlopesbyEvType, 'UniformOutput', false);
    findEmpt = cellfun(@isempty, semSlope);
    semSlope(findEmpt) = {0};
    semSlope = cell2mat(semSlope);
end %some empty bins

for sInd = 1:2
    subplot(1,2,sInd)
    lh = zeros(1,2);
    for g = 1:2
        hold on;
        lh(g) = errorbar(1:5, meanSlope(g,:,sInd), semSlope(g,:,sInd), semSlope(g,:,sInd), 'Color', rgb(cols{g}));
    end %group
    
    xlim([0 6])
    xticks(1:5)
    xlabel('Sleep session')
    
    ylim([100 400])
    ylabel('Normalize slope (°/event)')
    
    title([slopeNames{sInd} ' event slope (r^2 > ' num2str(r2Thresh) ')'])
    legend(lh, {group(1).name, group(2).name}, 'Location', 'northeastoutside')
end %ev direction

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 7 - PATH DIST

figtitle = 'ReplayPathDistance_acrossSleeps';
figure('Position', [717 504 525 388], 'Name', figtitle)

meanDist = cellfun(@mean, pathDist);
try
    semDist = cellfun(@semfunct, pathDist);
catch
    semDist = cellfun(@semfunct, pathDist, 'UniformOutput', false);
    findEmpt = cellfun(@isempty, semDist);
    semDist(findEmpt) = {0};
    semDist = cell2mat(semDist);
end %some empty bins

lh = zeros(1,2);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:5, meanDist(g,:), semDist(g,:), semDist(g,:), 'Color', rgb(cols{g}));
end %group

xlim([0 6])
xticks(1:5)
xlabel('Sleep session')

ylim([75 200])
ylabel('Path distance (°)')

title(['Replay event path distance (r^2 > ' num2str(r2Thresh) ')'])
legend(lh, {group(1).name, group(2).name})

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 8 - PATH DIST BY EV TYPE

figtitle = 'ReplayPathDistance_acrossSleeps_byEvType';
figure('Position', [472 504 1011 388], 'Name', figtitle)

meanDist = cellfun(@mean, pathDistbyEvType);
try
    semDist = cellfun(@semfunct, pathDistbyEvType);
catch
    semDist = cellfun(@semfunct, pathDistbyEvType, 'UniformOutput', false);
    findEmpt = cellfun(@isempty, semDist);
    semDist(findEmpt) = {0};
    semDist = cell2mat(semDist);
end %some empty bins

for sInd = 1:2
    subplot(1,2,sInd)
    lh = zeros(1,2);
    for g = 1:2
        hold on;
        lh(g) = errorbar(1:5, meanDist(g,:,sInd), semDist(g,:,sInd), semDist(g,:,sInd), 'Color', rgb(cols{g}));
    end %group
    
    xlim([0 6])
    xticks(1:5)
    xlabel('Sleep session')
    
    ylim([40 240])
    ylabel('Path distance (°)')
    
    title([slopeNames{sInd} ' event path distance (r^2 > ' num2str(r2Thresh) ')'])
    legend(lh, {group(1).name, group(2).name}, 'Location', 'northeastoutside')
end %ev direction

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end


%% FIG 9 - DURATION

figtitle = 'ReplayEventDuration_acrossSleeps';
figure('Position', [717 504 525 388], 'Name', figtitle)

meanDur = cellfun(@mean, eventDurs);
try
    semDur = cellfun(@semfunct, eventDurs);
catch
    semDur = cellfun(@semfunct, eventDurs, 'UniformOutput', false);
    findEmpt = cellfun(@isempty, semDur);
    semDur(findEmpt) = {0};
    semDur = cell2mat(semDur);
end %some empty bins

lh = zeros(1,2);
for g = 1:2
    hold on;
    lh(g) = errorbar(1:5, meanDur(g,:), semDur(g,:), semDur(g,:), 'Color', rgb(cols{g}));
end %group

xlim([0 6])
xticks(1:5)
xlabel('Sleep session')

ylim([0 0.2])
ylabel('Duration (s)')

title('All replay event duration')
legend(lh, {group(1).name, group(2).name})

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 8 - DURATION BY EV TYPE

figtitle = 'ReplayEventDuration_acrossSleeps_byEvType';
figure('Position', [472 504 1011 388], 'Name', figtitle)

meanDur = cellfun(@mean, evDursbyEvType);
try
    semDur = cellfun(@semfunct, evDursbyEvType);
catch
    semDur = cellfun(@semfunct, evDursbyEvType, 'UniformOutput', false);
    findEmpt = cellfun(@isempty, semDur);
    semDur(findEmpt) = {0};
    semDur = cell2mat(semDur);
end %some empty bins

for sInd = 1:2
    subplot(1,2,sInd)
    lh = zeros(1,2);
    for g = 1:2
        hold on;
        lh(g) = errorbar(1:5, meanDur(g,:,sInd), semDur(g,:,sInd), semDur(g,:,sInd), 'Color', rgb(cols{g}));
    end %group
    
    xlim([0 6])
    xticks(1:5)
    xlabel('Sleep session')
    
    ylim([0.05 0.25])
    ylabel('Duration (s)')
    
    title([slopeNames{sInd} ' event duration (r^2 > ' num2str(r2Thresh) ')'])
    legend(lh, {group(1).name, group(2).name}, 'Location', 'northeastoutside')
end %ev direction

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end



end %function