function fmr1CircTrack_x_lapByLap_runSpeed(group)
% function fmr1CircTrack_x_lapByLap_runSpeed(group)
%
% PURPOSE:
%   Plot run speed across laps for both genotypes.
%
% INPUT:
%   group struct
%
% OUTPUT:
%   Figure
%
% MMD
% 8/2021
% Colgin Lab

%% OPTIONS

saveOrNot = 0;
saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\BEHAVIOR';

prepForStats = 0;

%% INITIALIZE

minLapDay = inf;
for g = 1:2
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            dLaps = 0;
            for b = 1:4
                dLaps = dLaps + size(group(g).rat(r).day(d).begin(b).lapTms,1);
            end %begin
            if dLaps < minLapDay
                minLapDay = dLaps;
            end
        end %day
    end %rat
end %group

maxLapDay = 0;
for g = 1:2
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            dLaps = 0;
            for b = 1:4
                dLaps = dLaps + size(group(g).rat(r).day(d).begin(b).lapTms,1);
            end %begin
            if dLaps > maxLapDay
                maxLapDay = dLaps;
            end
        end %day
    end %rat
end %group

lapRunSpeed = cell(g,maxLapDay);

if prepForStats == 1
    statRunSpeed = [];
end

cols = {'Blue', 'Red'};

%% GET DATA

for g = 1:2
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            lapCntr = 0;
            
            if prepForStats == 1
                dayStat = [g zeros(1,minLapDay)];
            end %prep for stats
            for b = 1:4
                lapTms = group(g).rat(r).day(d).begin(b).lapTms;
                
                runSpeed = get_runspeed(group(g).rat(r).day(d).begin(b).coords);
                smSpeed = smooth_runspeed(runSpeed);
                for lp = 1:size(lapTms,1)
                    lapCntr = lapCntr + 1;
                    startInd = find(smSpeed(:,1) >= lapTms(lp,1), 1, 'First');
                    endInd = find(smSpeed(:,1) >= lapTms(lp,2), 1, 'First');
                    
                    try
                    lapRunSpeed{g,lapCntr} = [lapRunSpeed{g,lapCntr} mean(smSpeed(startInd:endInd,2))];
                    catch
                        keyboard
                    end
                    
                    if prepForStats == 1 && lapCntr <= minLapDay
                        dayStat(lapCntr+1) = mean(smSpeed(startInd:endInd,2));
                    end %stats
                end %lap
            end %begin
            if prepForStats == 1
                statRunSpeed = [statRunSpeed; dayStat];
            end %stats
        end %day
    end %rat
end %group
keyboard
%% PLOT FIGS

figtitle = 'RunSpeed_acrossLaps';

figure('Name', figtitle, 'Position', [594 467 699 420])

meanSpeed = cellfun(@mean, lapRunSpeed);
semSpeed = cellfun(@semfunct, lapRunSpeed, 'UniformOutput', false);

lh = nan(1,2);
leg = cell(2,1);
for g = 1:2
    % h = dotplot(1:minLapDay, lapRunSpeed(g,:), jitter, repmat(rgb(cols{g}),minLapDay,1), repmat(rgb('Black'),minLapDay,1));
    
   findSem = find(~cellfun(@isempty, semSpeed(g,:)));
   findMean = find(~isnan(meanSpeed(g,:)));
    
   if length(findSem) == length(findMean)
       lh(g) = errorbar(1:length(findMean), meanSpeed(g,1:findMean(end)), cell2mat(semSpeed(g,1:findMean(end))), cell2mat(semSpeed(g,1:findMean(end))), 'Color', rgb(cols{g})); 
   else
       semSpeed(g,findSem(end)+1:end) = {0};
       lh(g) = errorbar(1:length(findMean), meanSpeed(g,1:findMean(end)), cell2mat(semSpeed(g,1:findMean(end))), cell2mat(semSpeed(g,1:findMean(end))), 'Color', rgb(cols{g}));
   end
%     leg{g} = [group(g).name ' n = ' num2str(length(lapRunSpeed{g})) ' days'];
    hold on;
end %group

xlim([0 maxLapDay+1])
xticks(0:5:maxLapDay)
xlabel('Laps')

ylabel('Run speed (cm/s)')
% legend(lh, leg, 'Location', 'northeastoutside')
legend(group(1).name, group(2).name)

if saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end %saveOrNot

%% STATS?

if prepForStats == 1
    keyboard
end

end %function