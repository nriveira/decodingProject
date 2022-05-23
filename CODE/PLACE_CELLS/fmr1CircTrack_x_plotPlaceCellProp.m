function fmr1CircTrack_x_plotPlaceCellProp(group)
% function fmr1CircTrack_x_plotPlaceCellProp(group)
%
% PURPOSE:
%   Plot place cell firing properties as in Figure 3 of Mably et al. 2016
%   for data collected for the WT/FXS circle track data.
%
% INPUT:
%   group data struct
%
% OUTPUT:
%   Figures.
%
% OPTIONS:
%
% MM Donahue
% 5/2021
% Colgin Lab

%% OPTIONS

saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\PLACE_CELLS\placeCellProperties';
saveOrNot = 1;

prepForStats = 1;

cols = {'Blue', 'Red'};
groupNames = {'WT', 'FXS'};

spatBinSz = 4;

%% INITIALIZE

% Include all CA1 units, as in Mably et al. 2016 (Fig 3)
spatCorr = cell(2,4); %group x combo
rateOverlap = cell(2,4); %same

combos = [1 2; 2 3; 3 4; 1 4]; %combos as listed above
combosTxt = {'1-2', '2-3', '3-4', '1-4'};

% Cells with identified place fields
spatInfo = cell(2,1); %by group
meanFirRate = cell(2,1);
peakFieldFirRate = cell(2,1);
meanFieldFirRate = cell(2,1);
fieldSize = cell(2,1);
fieldPerCell = cell(2,1);

convFact = (2*pi*50)/(360/spatBinSz); %cm per bin, because track is 100 cm in diameter and ratemap has 72 bins

uCntr = 0; %initialize for stats (ANOVAs)


%% GET DATA

for g = 1:2
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            for u = 1:length(group(g).rat(r).day(d).xBeginUnitInfo)
                uCntr = uCntr + 1;
                uRateMaps = zeros(4,360/spatBinSz); %for storing rate maps from all begins for this unit
                
                for b = 1:4
                    uRateMaps(b,:) = group(g).rat(r).day(d).begin(b).unit(u).smRateMap; %get all of the ratemaps for the spat corr comparisons
                end %begin
                
                for c = 1:4
                    b1 = combos(c,1); %begins to compare
                    b2 = combos(c,2);
                    
                    rm1 = uRateMaps(b1,:); %ratemaps to compare
                    rm2 = uRateMaps(b2,:);
                    
                    scMat = corrcoef(rm1,rm2);
                    if ~isnan(scMat(2))
                        spatCorr{g,c} = [spatCorr{g,c} scMat(2)];
                    end

                    fr1 = mean(uRateMaps(b1,:));
                    fr2 = mean(uRateMaps(b2,:));
                    
                    rateRatio = min([fr1 fr2]) / max([fr1 fr2]); % Calculate the rate ratio with the lower FR as numerator
                    if ~isnan(rateRatio)
                        rateOverlap{g,c} = [rateOverlap{g,c} rateRatio];
                    end
                    
                end %combos
                
                if ~isempty(group(g).rat(r).day(d).xBeginUnitInfo(u).pf) %only cells with place fields from avg ratemap
                    
                    rateMap = group(g).rat(r).day(d).xBeginUnitInfo(u).rateMap; %get AVG ratemap across all begins
                    timePerBin = group(g).rat(r).day(d).xBeginTpb; %same, all begins
                    
                    tmpSpatInfo = get_spatial_info(rateMap, timePerBin);
                    spatInfo{g,1} = [spatInfo{g,1} tmpSpatInfo];
                    
                    fieldPerCell{g,1} = [fieldPerCell{g,1} length(group(g).rat(r).day(d).xBeginUnitInfo(u).pf)];
                    
                    tmpPeak = [];
                    tmpFieldMean = [];
                    tmpSize = [];
                    
                    for pf = 1:length(group(g).rat(r).day(d).xBeginUnitInfo(u).pf) %if there are multiple place fields, get mean of measures for this unit
                        tmpPeak = [tmpPeak group(g).rat(r).day(d).xBeginUnitInfo(u).pf(pf).pkFr];
                        
                        pfInds = group(g).rat(r).day(d).xBeginUnitInfo(u).pf(pf).inds;
                        pullBins = zeros(1,length(pfInds)); %initialize so I know it works right
                        
                        if pfInds(1) < pfInds(end) %it doesn't wrap around 0 degrees
                            pullBins = rateMap(1,pfInds(1):pfInds(end));
                        else %pf does wrap around 0 degrees
                            findWrap = find(pfInds == 1); %find where it crosses 0 degree boundary
                            pullBins(1,1:findWrap-1) = rateMap(1,pfInds(1):pfInds(findWrap-1));
                            pullBins(1,findWrap:end) = rateMap(1,pfInds(findWrap):pfInds(end));
                        end %whether or not place field wraps around 0
                        
                        tmpFieldMean = [tmpFieldMean mean(pullBins)];
                        
                        tmpSize = [tmpSize length(pfInds)*convFact];
                        
                    end %pf
                    
                    peakFieldFirRate{g,1} = [peakFieldFirRate{g,1} mean(tmpPeak)];
                    meanFieldFirRate{g,1} = [meanFieldFirRate{g,1} mean(tmpFieldMean)];
                    
                    meanFirRate{g,1} = [meanFirRate{g,1} mean(uRateMaps(:))];
                    fieldSize{g,1} = [fieldSize{g,1} mean(tmpSize)];
                    
                end %if this unit has a place field
                
            end %unit
        end %day
    end %rat
    
end %group

%% MAKE FIGURES

% SPATIAL CORRELATION
keyboard
avgSpatCorr = cellfun(@mean, spatCorr);
semSpatCorr = cellfun(@semfunct, spatCorr);

figtitle = ['SpatialCorrelation_allCells'];

figure('Name', figtitle)

lh = NaN(1,2);
for g = 1:2
    hold on;
    lh(g) =  errorbar(1:4, avgSpatCorr(g,:), semSpatCorr(g,:), 'Color', cols{g}, 'LineWidth', 1); %just get the line handle for the first one for the legend
    %     errorbar(4:5, avgSpatCorr(g,4:5), semSpatCorr(g,4:5), 'Color', cols{g}, 'LineWidth', 1);
    %     errorbar(6, avgSpatCorr(g,6), semSpatCorr(g,6), 'Color', cols{g}, 'LineWidth', 1);
    
end %group

ylabel('Spatial correlation')
ylim([0 1])
xlim([0 5])
xticks(1:4)
xticklabels(combosTxt)
xlabel('Sessions compared')
legend(lh, groupNames, 'Location', 'southeast')



if saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

% RATE OVERLAP

avgRateOverlap = cellfun(@mean, rateOverlap);
semRateOverlap = cellfun(@semfunct, rateOverlap);

figtitle = ['RateOverlap_allCells'];

figure('Name', figtitle)

lh = NaN(1,2);
for g = 1:2
    hold on;
    lh(g) =  errorbar(1:4, avgRateOverlap(g,:), semRateOverlap(g,:), 'Color', cols{g}, 'LineWidth', 1); %just get the line handle for the first one for the legend
    %     errorbar(4:5, avgRateOverlap(g,4:5), semRateOverlap(g,4:5), 'Color', cols{g}, 'LineWidth', 1);
    %     errorbar(6, avgRateOverlap(g,6), semRateOverlap(g,6), 'Color', cols{g}, 'LineWidth', 1);
    
end %group

ylabel('Rate overlap')
ylim([0 1])
xlim([0 5])
xticks(1:4)
xticklabels(combosTxt)
xlabel('Sessions compared')
legend(lh, groupNames, 'Location', 'southeast')

% fprintf('Rate overlap mixed repeated measures ANOVA:\n')
% mixed_between_within_anova(statRateOverlap);

if saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

% PLACE CELL INFO

figtitle = ['PlaceCellInfo_allCells'];

figure('Name', figtitle, 'Position', [146 274 1269 674])

% Spatial information
subplot(2,3,1)

avgSpatInfo = cellfun(@mean, spatInfo);

xVals = 1:2;
jitter = 0.01;
h = dotplot(1:2, spatInfo, jitter, [rgb(cols{1}); rgb(cols{2})], [0 0 0; 0 0 0]);

for g = 1:2
    line([g-.3 g+.3], [avgSpatInfo(g) avgSpatInfo(g)], 'Color', 'black', 'LineWidth', 3)
end %group

ylabel('Spatial information (bits/spike)')
ylim([0 10])
xticklabels(groupNames)

[P,~,STATS] = ranksum(spatInfo{1,1}, spatInfo{2,1});
fprintf('Spatial information:\n')
fprintf('\tMann-Whitney U = %.2g, n1 = %d, n2 = %d, p = %.3g\n', STATS.ranksum, length(spatInfo{1,1}), length(spatInfo{2,1}), P);
fprintf('\t\tMean - WT: %.2g; KO: %.2g\n', mean(spatInfo{1,1}), mean(spatInfo{2,1}));

if P < 0.05
    yLimits = ylim;
    line([1 2], [yLimits(2) - 1 yLimits(2) - 1], 'Color', 'black', 'LineWidth', 1);
    text(1.25, yLimits(2) - 0.5, ['p = ' num2str(round(P,3))])
end %significant p value

% Place field size

subplot(2,3,2)

avgFieldSize = cellfun(@mean, fieldSize);

h = dotplot(1:2, fieldSize, jitter, [rgb(cols{1}); rgb(cols{2})], [0 0 0; 0 0 0]);

for g = 1:2
    line([g-.3 g+.3], [avgFieldSize(g) avgFieldSize(g)], 'Color', 'black', 'LineWidth', 3)
end %group

ylabel('Place field size (cm)')
ylim([0 120])
xticklabels(groupNames)

[P,~,STATS] = ranksum(fieldSize{1,1}, fieldSize{2,1});
fprintf('Place field size:\n')
fprintf('\tMann-Whitney U = %.2g, n1 = %d, n2 = %d, p = %.3g\n', STATS.ranksum, length(fieldSize{1,1}), length(fieldSize{2,1}), P);
fprintf('\t\tMedians - WT: %.2g; KO: %.2g\n', median(fieldSize{1,1}), median(fieldSize{2,1}));

if P < 0.05
    yLimits = ylim;
    line([1 2], [yLimits(2) - 15 yLimits(2) - 15], 'Color', 'black', 'LineWidth', 1);
    text(1.25, yLimits(2) - 10, ['p = ' num2str(round(P,3))])
end %significant p value

% Place fields per cell

subplot(2,3,3)

avgFieldPerCell = cellfun(@mean, fieldPerCell);
semFieldPerCell = cellfun(@semfunct, fieldPerCell);

bgraph = bar(avgFieldPerCell);
bgraph.FaceColor = 'flat';
for g = 1:2 %both types
    bgraph.CData(g,:) = rgb(cols{g});
end %types

hold on;
er = errorbar(1:2, avgFieldPerCell, semFieldPerCell, semFieldPerCell);
er.Color = [0 0 0];
er.LineStyle = 'none';

ylabel('Place fields/cell')
ylim([0 1.5])
xticklabels(groupNames)

[P,~,STATS] = ranksum(fieldPerCell{1,1}, fieldPerCell{2,1});
fprintf('Place fields per cell:\n')
fprintf('\tMann-Whitney U = %.2g, n1 = %d, n2 = %d, p = %.3g\n', STATS.ranksum, length(fieldPerCell{1,1}), length(fieldPerCell{2,1}), P);

if P < 0.05
    yLimits = ylim;
    line([1 2], [yLimits(2) - .15 yLimits(2) - .15], 'Color', 'black', 'LineWidth', 1);
    text(1.25, yLimits(2) - 0.1, ['p = ' num2str(round(P,3))])
end %significant p value


%  Mean firing rate

subplot(2,3,4)

avgMeanFirRate = cellfun(@mean,meanFirRate);

h = dotplot(1:2, meanFirRate', jitter, [rgb(cols{1}); rgb(cols{2})], [0 0 0; 0 0 0]);

for g = 1:2
    line([g-.3 g+.3], [avgMeanFirRate(g) avgMeanFirRate(g)], 'Color', 'black', 'LineWidth', 3)
end %group

ylabel('Mean firing rate (Hz)')
xticklabels(groupNames)

[P,~,STATS] = ranksum(meanFirRate{1,1}, meanFirRate{2,1});
fprintf('Mean firing rate:\n')
fprintf('\tMann-Whitney U = %.2g, n1 = %d, n2 = %d, p = %.3g\n', STATS.ranksum, length(meanFirRate{1,1}), length(meanFirRate{2,1}), P);
fprintf('\t\tMedians - WT: %.2g; KO: %.2g\n', median(meanFirRate{1,1}), median(meanFirRate{2,1}));

if P < 0.05
    yLimits = ylim;
    line([1 2], [yLimits(2) - 1 yLimits(2) - 1], 'Color', 'black', 'LineWidth', 1);
    text(1.25, yLimits(2) - 0.5, ['p = ' num2str(round(P,3))])
end %significant p value

%Peak firing rate

subplot(2,3,5)

avgPeakFirRate = cellfun(@mean, peakFieldFirRate);
h = dotplot(1:2, peakFieldFirRate, jitter, [rgb(cols{1}); rgb(cols{2})], [0 0 0; 0 0 0]);

for g = 1:2
    line([g-.3 g+.3], [avgPeakFirRate(g) avgPeakFirRate(g)], 'Color', 'black', 'LineWidth', 3)
end %group

ylabel('Peak in-field firing rate (Hz)')
xticklabels(groupNames)

[P,~,STATS] = ranksum(peakFieldFirRate{1,1}, peakFieldFirRate{2,1});
fprintf('Peak in-field firing rate:\n')
fprintf('\tMann-Whitney U = %.2g, n1 = %d, n2 = %d, p = %.3g\n', STATS.ranksum, length(peakFieldFirRate{1,1}), length(peakFieldFirRate{2,1}), P);
fprintf('\t\tMedians - WT: %.2g; KO: %.2g\n', median(peakFieldFirRate{1,1}), median(peakFieldFirRate{2,1}));

if P < 0.05
    yLimits = ylim;
    line([1 2], [yLimits(2) - 5 yLimits(2) - 5], 'Color', 'black', 'LineWidth', 1);
    text(1.25, yLimits(2) - 3, ['p = ' num2str(round(P,3))])
end %significant p value

% Average in-field firing rate

subplot(2,3,6)

avgMeanFieldFirRate = cellfun(@mean, meanFieldFirRate);
h = dotplot(1:2, meanFieldFirRate, jitter, [rgb(cols{1}); rgb(cols{2})], [0 0 0; 0 0 0]);

for g = 1:2
    line([g-.3 g+.3], [avgMeanFieldFirRate(g) avgMeanFieldFirRate(g)], 'Color', 'black', 'LineWidth', 3)
end %group

ylabel('Mean in-field firing rate (Hz)')
xticklabels(groupNames)

[P,~,STATS] = ranksum(meanFieldFirRate{1,1}, meanFieldFirRate{2,1});
fprintf('Mean in-field firing rate:\n')
fprintf('\tMann-Whitney U = %.2g, n1 = %d, n2 = %d, p = %.3g\n', STATS.ranksum, length(meanFieldFirRate{1,1}), length(meanFieldFirRate{2,1}), P);
fprintf('\t\tMedians - WT: %.2g; KO: %.2g\n', median(meanFieldFirRate{1,1}), median(meanFieldFirRate{2,1}));


if P < 0.05
    yLimits = ylim;
    line([1 2], [yLimits(2) - 5 yLimits(2) - 5], 'Color', 'black', 'LineWidth', 1);
    text(1.25, yLimits(2) - 3, ['p = ' num2str(round(P,3))])
end %significant p value

if saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

%% STATS?

if prepForStats == 1
    
    statSpatCorr = [];
    statRateOverlap = [];
    statSpatInfo = [];
    statFieldSize = [];
    statMeanFirRate = [];
    statPeakFieldFirRate = [];
    statMeanFieldFirRate = [];
    statFieldPerCell = [];
    
    for g = 1:2
        c = 1; %temporarily
        numU = length(spatCorr{g,c});
        for u = 1:numU
            tmpStatSC = g;
            tmpStatRO = g;
            
            for c = 1:4
                tmpStatSC = [tmpStatSC spatCorr{g,c}(u)];
                tmpStatRO = [tmpStatRO rateOverlap{g,c}(u)];
            end %combos
            statSpatCorr = [statSpatCorr; tmpStatSC];
            statRateOverlap = [statRateOverlap; tmpStatRO];
            
        end %unit
        
        for u = 1:length(spatInfo{g})
            statSpatInfo = [statSpatInfo; g spatInfo{g}(u)];
            
            statFieldSize = [statFieldSize; g fieldSize{g}(u)];
            
            statMeanFirRate = [statMeanFirRate; g meanFirRate{g}(u)];
            statPeakFieldFirRate = [statPeakFieldFirRate; g peakFieldFirRate{g}(u)];
            statMeanFieldFirRate = [statMeanFieldFirRate; g meanFieldFirRate{g}(u)];
            
            statFieldPerCell = [statFieldPerCell fieldPerCell{g}(u)];
        end %u - spat info
        
    end %group
    
    keyboard
end %prep for stats





end %function


