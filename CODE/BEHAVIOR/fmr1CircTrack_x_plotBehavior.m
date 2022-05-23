function fmr1CircTrack_x_plotBehavior(group)
% function fmr1CircTrack_x_plotBehavior(group)
%
% PURPOSE:
%   Examine any potenital differences in behavior between WT and KO groups.
%
% INPUT:
%   group struct
%
% OUTPUT:
%   Figures:
%       F1: Run speed means overall and across begins.
%       F2: Laps completed overall and across beegins.
%
% MMD
% 7/2021
% Colgin Lab

%% OPTIONS

saveOrNot = 1;

%% INITIALIZE

saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\behavior';
curDir = pwd;

runSpeedAll = cell(2,1); %by group
runSpeedByBegin = cell(2,4); %by group and begin

numLapsAll = cell(2,1);
numLapsByBegin = cell(2,4);

cols = {'Blue', 'Red'};

%% GET DATA

for g = 1:2
    for r = 1:length(group(g).rat)
        %         ratSpeedAll = [];
        %         ratSpeedByBegin = cell(1,4);
        %
        %         ratLapsAll = [];
        %         ratLapByBegin = cell(1,4);
        
        for d = 1:length(group(g).rat(r).day)
            daySpeed = [];
            dayLaps = [];
            for b = 1:4
                
                coords = group(g).rat(r).day(d).begin(b).coords;
                instRs = get_runspeed(coords);
                smRs = smooth_runspeed(instRs);
                if isnan(mean(smRs))
                    keyboard
                end
                %                 ratSpeedAll = [ratSpeedAll mean(smRs(:,2))];
                %                 ratSpeedByBegin{b} = [ratSpeedByBegin{b} mean(smRs(:,2))];
                daySpeed = [daySpeed mean(smRs(:,2))];
                runSpeedByBegin{g,b} = [runSpeedByBegin{g,b} mean(smRs(:,2))];
                
                numLaps = size(group(g).rat(r).day(d).begin(b).lapInds,1);
                %                 ratLapsAll = [ratLapsAll numLaps];
                dayLaps = [dayLaps numLaps];
                numLapsByBegin{g,b} = [numLapsByBegin{g,b} numLaps];
                %                 ratLapByBegin{b} = [ratLapByBegin{b} numLaps];
                
            end %begin
            runSpeedAll{g} = [runSpeedAll{g} mean(daySpeed)];
            numLapsAll{g} = [numLapsAll{g} mean(dayLaps)];
            
        end %day
        
        %         meanSpeedByBegin = cellfun(@mean, ratSpeedByBegin);
        %         meanLapsByBegin = cellfun(@mean, ratLapByBegin);
        %         runSpeedAll{g} = [runSpeedAll{g} mean(ratSpeedAll)];
        %         numLapsAll{g} = [numLapsAll{g} mean(ratLapsAll)];
        
        %         for b = 1:4
        %             runSpeedByBegin{g,b} = [runSpeedByBegin{g,b} meanSpeedByBegin(b)];
        %             numLapsByBegin{g,b} = [numLapsByBegin{g,b} meanLapsByBegin(b)];
        %         end
        
    end %rat
end %group


%% FIG 1

figtitle = 'RunSpeed';

figure('Name', figtitle, 'Position', [461 456 1039 420])

subplot(1,2,1)

jitter = 0.1;
dotplot(1:2, runSpeedAll, jitter, [rgb(cols{1}); rgb(cols{2})], [0 0 0; 0 0 0]);
ylabel('Run speed (cm/s)')


meanSpd = cellfun(@mean, runSpeedAll);
xLabs = cell(1,2);
for g = 1:2
    line([g-0.15 g+0.15], [meanSpd(g) meanSpd(g)], 'Color', 'Black', 'LineWidth', 3)
    
    xLabs{g} = [group(g).name ' n = ' num2str(length(runSpeedAll{g})) ' days'];
end

xticklabels(xLabs)
ylim([30 55])

subplot(1,2,2)
jitter = 0.05; %smaller for this more confusing plot
for g = 1:2
    faceColors = repmat(rgb(cols{g}),4,1); %jsut make it easier to plug into the formula
    edgeColors = repmat(rgb('Black'),4,1);
    
    h = dotplot(1:4, runSpeedByBegin(g,:), jitter, faceColors, edgeColors);
    hold on;
end

meanSpdBeg = cellfun(@mean, runSpeedByBegin);
lh = nan(1,2);
for g = 1:2
    for b = 1:4
        line([b-0.15 b+0.15], [meanSpdBeg(g,b) meanSpdBeg(g,b)], 'Color', rgb(cols{g}), 'LineWidth', 3)
        
    end %begin
    
    lh(g) = plot(1:4, meanSpdBeg(g,:), 'Color', rgb(cols{g}));
end %group

ylabel('Run speed (cm/s)')
ylim([30 55])
xlabel('Begin #')
legend(lh, xLabs, 'Location', 'southeast')

if  saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    cd(curDir)
end

%% FIG 2 LAPS

figtitle = 'NumberofLaps';

figure('Name', figtitle, 'Position', [461 456 1039 420])

subplot(1,2,1)

jitter = 0.1;
dotplot(1:2, numLapsAll, jitter, [rgb(cols{1}); rgb(cols{2})], [0 0 0; 0 0 0]);
ylabel('Number of laps/ession')

meanLaps = cellfun(@mean, numLapsAll);
xLabs = cell(1,2);
for g = 1:2
    line([g-0.15 g+0.15], [meanLaps(g) meanLaps(g)], 'Color', 'Black', 'LineWidth', 3)
    
    xLabs{g} = [group(g).name ' n = ' num2str(length(numLapsAll{g})) ' days'];
end

xticklabels(xLabs)
ylim([3 12])

subplot(1,2,2)
jitter = 0.05; %smaller for this more confusing plot
for g = 1:2
    faceColors = repmat(rgb(cols{g}),4,1); %jsut make it easier to plug into the formula
    edgeColors = repmat(rgb('Black'),4,1);
    
    h = dotplot(1:4, numLapsByBegin(g,:), jitter, faceColors, edgeColors);
    hold on;
end

meanLapBeg = cellfun(@mean, numLapsByBegin);
lh = nan(1,2);
for g = 1:2
    for b = 1:4
        line([b-0.15 b+0.15], [meanLapBeg(g,b) meanLapBeg(g,b)], 'Color', rgb(cols{g}), 'LineWidth', 3)
        
    end %begin
    
    lh(g) = plot(1:4, meanLapBeg(g,:), 'Color', rgb(cols{g}));
end %group

ylabel('Laps/session')
ylim([0 12])
xlabel('Begin #')
legend(lh, xLabs, 'Location', 'southeast')

if  saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    cd(curDir)
end

end %function