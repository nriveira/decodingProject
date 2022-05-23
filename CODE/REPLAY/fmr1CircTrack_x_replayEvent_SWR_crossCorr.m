function fmr1CircTrack_x_replayEvent_SWR_crossCorr(group)
% function fmr1CircTrack_x_replayEvent_SWR_crossCorr(group)
% 
% PURPOSE:
%   Compare the detected replay event (from cell population activity) times
%   to the detected SWR (from the LFP) times.
% 
% INPUT:
%   group = data struct, through fxn 6
% 
% OUTPUT:
%   Figures showing cross correlation of replay event and detected SWR
%   times. See fmr1CircTrack_6_addReplayEventsToStruct for info on how SWR
%   times are adjusted.
%       F1: Non-adjusted SWR times
%       F2: Adjusted SWR times
% 
% MMD
% 11/2021
% Colgin Lab

%% OPTIONS

saveOrNot = 1;

binSz = 0.1; %.1 s
timeLag = 0.5;

%% INITIALIZE

nonAdjxCorr = cell(2,1); %by group - non-adjusted SWR times
adjxCorr = cell(2,1); %adjusted SWR times

cols = {'Blue', 'Red'};

saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\SHARP_WAVE_RIPPLES\CROSS_CORR';
curDir = pwd;
cd(saveDir);

%% GET DATA

for g = 1:2
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            for s = 1:5
                if isempty(group(g).rat(r).day(d).sleep(s).unit) %if there is anything in the sleep
                    continue % to next sleep
                end %anything in sleep
                
                replayEvTms = zeros(1, length(group(g).rat(r).day(d).sleep(s).popEv));
                for i = 1:length(group(g).rat(r).day(d).sleep(s).popEv)
                    replayEvTms(i) = group(g).rat(r).day(d).sleep(s).popEv(i).tms(1);
                end %i - pop events
                
                SWRtms = zeros(1, length(group(g).rat(r).day(d).sleep(s).rip));
                adjSWRtms = nan(1, length(group(g).rat(r).day(d).sleep(s).rip));
                for i = 1:length(group(g).rat(r).day(d).sleep(s).rip)
                    SWRtms(i) = group(g).rat(r).day(d).sleep(s).rip(i).tms(1);
                    if ~isempty(group(g).rat(r).day(d).sleep(s).rip(i).adjTms)
                        adjSWRtms(i) = group(g).rat(r).day(d).sleep(s).rip(i).adjTms(1);
                    end %there are adjusted times
                end %i - rip
                
                slpStart = group(g).rat(r).day(d).sleep(s).coords(1,1);
                slpEnd = group(g).rat(r).day(d).sleep(s).coords(end,1);
                
                edges = slpStart:binSz:slpEnd; %get bin edges
                replayBinCnts = histcounts(replayEvTms, edges);
                SWRbinCnts = histcounts(SWRtms, edges);
                adjSWRbinCnts = histcounts(adjSWRtms, edges);
                
                if ~isempty(replayEvTms) && ~isempty(SWRtms) && ~isempty(adjSWRtms)
                    tmpCorr = xcorr(replayBinCnts, SWRbinCnts,  timeLag/binSz, 'coeff');
                    nonAdjxCorr{g} = [nonAdjxCorr{g}; tmpCorr];
                    
                    tmpAdjCorr = xcorr(replayBinCnts, adjSWRbinCnts,  timeLag/binSz, 'coeff');
                    adjxCorr{g} = [adjxCorr{g}; tmpAdjCorr];
                end %thre is anything to correlate
                
            end %sleep
        end %day
    end %rat
end %group


%% FIG 1 - NON AJUSTED TIMES

figtitle = 'ReplayEvent_SWR_nonAdj_xCorr';
figure('Name', figtitle, 'Position', [559 545 731 420])

tmpMean = cellfun(@mean, nonAdjxCorr, 'UniformOutput', false);
tmpStd = cellfun(@semfunct, nonAdjxCorr, 'UniformOutput', false);

yMax = 0;
for g = 1:2
    subplot(1,2,g)
    bgraph = bar(-timeLag:binSz:timeLag, tmpMean{g}, 'FaceColor', 'Flat');
    bgraph.CData(:) = repmat(rgb(cols{g}), length(tmpMean{g}), 1);
    
    hold on;
    errorbar(-timeLag:binSz:timeLag, tmpMean{g}, tmpStd{g}, tmpStd{g}, 'Color', 'Black', 'LineStyle', 'None')
    
    xlabel('Time lag (s)')
    ylabel('Proportion of non-adjusted SWR events')
    title(group(g).name)
    
    ax = gca;
    if ax.YLim(2) > yMax
        yMax = ax.YLim(2);
    end %store yMax
end %group

for g = 1:2
    ylim([0 yMax])
end %group

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %save or not

%% FIG 2 - ADJUSTED TIMES

figtitle = 'ReplayEvent_SWR_adj_xCorr';
figure('Name', figtitle, 'Position', [559 545 731 420])

tmpMean = cellfun(@mean, adjxCorr, 'UniformOutput', false);
tmpStd = cellfun(@semfunct, adjxCorr, 'UniformOutput', false);

yMax = 0;
for g = 1:2
    subplot(1,2,g)
    bgraph = bar(-timeLag:binSz:timeLag, tmpMean{g}, 'FaceColor', 'Flat');
    bgraph.CData(:) = repmat(rgb(cols{g}), 11, 1);
    
    hold on;
    errorbar(-timeLag:binSz:timeLag, tmpMean{g}, tmpStd{g}, tmpStd{g}, 'Color', 'Black', 'LineStyle', 'None')
    
    xlabel('Time lag (s)')
    ylabel('Proportion of adjusted SWR events')
    title(group(g).name)
    
    ax = gca;
    if ax.YLim(2) > yMax
        yMax = ax.YLim(2);
    end %store yMax
end %group

for g = 1:2
    ylim([0 yMax])
end %group

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %save or not

cd(curDir)

end %function