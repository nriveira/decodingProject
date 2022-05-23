function fmr1CircTrack_x_replayEvent_autoCorr(group)
% function fmr1CircTrack_x_replayEvent_autoCorr(group)
%



%% OPTIONS

saveOrNot = 1;

binSz = .1;% in s
timeLag = 5;

%% INITIALIZE

evAutoCorr = cell(2,1); %by group

cols = {'Blue', 'Red'};

saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\REPLAY\replayRate';
curDir = pwd;
cd(saveDir);

%% GET DATA


for g = 1:2
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            for s = 1:5
                if isempty(group(g).rat(r).day(d).sleep(s).unit) || isempty(group(g).rat(r).day(d).sleep(s).popEv) %if there is anything in the sleep
                    continue % to next sleep - not break some just don't have sleep1
                end %anything in sleep
                
                replayEvTms = zeros(1, length(group(g).rat(r).day(d).sleep(s).popEv));
                for i = 1:length(group(g).rat(r).day(d).sleep(s).popEv)
                    replayEvTms(i) = group(g).rat(r).day(d).sleep(s).popEv(i).tms(1);
                end %i - pop events

                slpStart = group(g).rat(r).day(d).sleep(s).coords(1,1);
                slpEnd = group(g).rat(r).day(d).sleep(s).coords(end,1);
                
                edges = slpStart:binSz:slpEnd; %get bin edges
                
                replayBinCnts = histcounts(replayEvTms, edges);
                
                tmpCorr = xcorr(replayBinCnts,  timeLag/binSz, 'coeff');
                evAutoCorr{g} = [evAutoCorr{g}; tmpCorr];
                
            end %sleep
        end %day
    end %rat
end %group
keyboard
%% FIGS

figtitle = 'ReplayEvent_autoCorr';
figure('Name', figtitle, 'Position', [559 545 731 420])

tmpMean = cellfun(@mean, evAutoCorr, 'UniformOutput', false);
tmpStd = cellfun(@semfunct, evAutoCorr, 'UniformOutput', false);

yMax = 0;
for g = 1:2
    subplot(1,2,g)
    bgraph = bar(-timeLag:binSz:timeLag, tmpMean{g}, 'FaceColor', 'Flat');
    bgraph.CData(:) = repmat(rgb(cols{g}), length(tmpMean{g}), 1);
    
    hold on;
%     errorbar(-timeLag:binSz:timeLag, tmpMean{g}, tmpStd{g}, tmpStd{g}, 'Color', 'Black', 'LineStyle', 'None')
    
    xlabel('Time lag (s)')
    ylabel('Proportion replay events')
    title(group(g).name)
    
    ylim([0 0.1])

end %group

end %function