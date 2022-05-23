function fmr1CircTrack_x_shortTermStability(group)
% function fmr1CircTrack_x_shortTermStability(group)
% 
% PURPOSE:
%   Examine short-term stability of place cells by computing spatial
%   correlation between first and last half of each begin. Inspired by
%   differences seen in short-term stability in WT and KO mice in Arbab,
%   Pennartz et al. 2018. The 10-min begin sessions for circle track are
%   pretty short so I'm not sure how robust or useful this will be.
%
% INPUT:
%   group data struct
%
% OUTPUT:
%   Figure.
% 
% NOTE:
%   Complementary function fmr1CircTrack_x_lapRasterPlots(group) could be
%   used for illustrative purposes if there are differences in short term
%   stability between groups.
%
% MM Donahue
% 5/2021
% Colgin Lab


%% OPTIONS

saveOrNot = 1;
saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\PLACE_CELLS\placeCellProperties';

prepForStats = 1;

cols = {'Blue', 'Red'};
groupNames = {'WT', 'FXS'};

%% INITIALIZE

shortTermCorr = cell(2,4); %group x session

%for get_ratemap code
spatBinSz = 5; %degrees
velFilt = 1;
durCrit = 1;

jitter = 0.1; %for dotplot

curDir = pwd;
%% GET DATA

for g = 1:2
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            for b = 1:4
                radPos = group(g).rat(r).day(d).begin(b).radPos;
                coords = group(g).rat(r).day(d).begin(b).coords; %going to use it a lot
                
                begStart = coords(1,1);
                begEnd = coords(end,1);
                
                halfInd = round(size(coords,1)/2);
                begHalf = coords(halfInd,1);
                
                for u = 1:length(group(g).rat(r).day(d).xBeginUnitInfo)
                    
                    spkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                    
                    spkTmsFirst = find(spkTms < begHalf);
                    spkTmsLast = find(spkTms >= begHalf); 
                    
                    rateMapFirst = get_ratemap_circtrack(spkTms(spkTmsFirst), coords(1:halfInd,:), radPos(1:halfInd,:), spatBinSz, velFilt, durCrit);
                    rateMapLast = get_ratemap_circtrack(spkTms(spkTmsLast), coords(halfInd:end,:), radPos(halfInd:end,:), spatBinSz, velFilt, durCrit);
                    
                    scMat = corrcoef(rateMapFirst, rateMapLast);
                    
                    shortTermCorr{g,b} = [shortTermCorr{g,b} scMat(2)];
                    
                    
                end %unit
            end %begin

        end %day
    end %rat
end %group
keyboard

%% FIGS

figtitle = 'ShortTermStability';

figure('Name', figtitle, 'Position', [396 529 1158 449])
% xValOpts = [0.85:1:3.85; 1.15:1:4.15];
% 
% subplot(1,2,1)
% lh = zeros(1,2);
% legText = cell(2,1);
% for g = 1:2
%     xVals = xValOpts(g,:);
%     faceColors = repmat(rgb(cols{g}),4,1);
%     edgeColors = repmat(rgb('Black'),4,1);
%     h = dotplot(xVals, shortTermCorr(g,:), jitter, faceColors, edgeColors);
%     lh(g) = h(1);
%     
%     legText{g} = [groupNames{g} ': n = ' num2str(length(shortTermCorr{g,1})) ' cells'];
%     
%     tmpMean = cellfun(@nanmean, shortTermCorr(g,:));
%     for b = 1:4
%         line([xValOpts(g,b)-0.1 xValOpts(g,b)+0.1], [tmpMean(b) tmpMean(b)], 'Color', rgb('Black'), 'LineWidth', 3);
%     end %begin
% end %group
% 
% xlim([0 5])
% xticks(1:4)
% xlabel('Session')
% ylim([0 1])
% ylabel('Spatial correlation')
% legend(lh, legText, 'Location', 'southoutside')
% title('Spatial correlation of ratemaps from first and second half of session')

subplot(1,2,1)

tmpMean = cellfun(@nanmean, shortTermCorr);
tmpSEM = cellfun(@nansemfunct, shortTermCorr);

lh = NaN(1,2);
for g = 1:2
    hold on;
    
    lh(g) =  errorbar(1:4, tmpMean(g,:), tmpSEM(g,:), 'Color', cols{g}, 'LineWidth', 1);
    
     legText{g} = [groupNames{g} ': n = ' num2str(length(shortTermCorr{g,1})) ' cells'];
end %group

xlim([0 5])
xticks(1:4)
xlabel('Session')
ylim([0 0.7])
ylabel('Spatial correlation')
legend(lh, legText, 'Location', 'southwest')
title('Spatial correlation of ratemaps from first and second half of session')


subplot(1,2,2)

tmpMean = tmpMean';
tmpSEM = tmpSEM';

bgraph = bar(1:4, tmpMean, 'FaceColor', 'Flat');
hold on;

ngroups = 4;
nbars = 2;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for g = 1:2
    bgraph(g).CData = rgb(cols{g});
    
    x = (1:ngroups) - groupwidth/2 + (2*g-1) * groupwidth / (2*nbars);
    er = errorbar(x, tmpMean(:,g), tmpSEM(:,g), '.');
    %     er = errorbar(1:4, tmpMean(g,:), tmpSEM(g,:), tmpSEM(g,:));
    er.Color = [0 0 0];
    %     er.LineStyle = 'none';
end %group

xticks(1:4)
xlabel('Session')
ylabel('Spatial correlation')
ylim([0 0.7])
legend(legText, 'Location', 'northeastoutside')
% title('Spatial correlation of ratemaps from first and second half of session')


if saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
    cd(curDir)
end %save option

%% STATS

if prepForStats == 1
    
    statShortTermCorr = [];
    
    for g = 1:2
        for u = 1:length(shortTermCorr{g,1}) %same num of cells in all sessions for each group
            
            tmpStat = g;
            for b = 1:4
                tmpStat = [tmpStat shortTermCorr{g,b}(u)];
            end %begin
            statShortTermCorr = [statShortTermCorr; tmpStat];
            
        end %unit
    end %group
    
    keyboard
end %stats


end %function