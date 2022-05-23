function fmr1CircTrack_x_phase_precession(group)
% function fmr1CircTrack_x_phase_precession(group)
%
% PURPOSE:
%   Plot phase precession of neurons on the circle track from WT and KO
%   groups.
%
% INPUT:
%   group struct
%
% OUTPUT:
%   Figures:
%       F1: Heat plots of spike counts for theta phase x normalized
%           position in place field (as in Bieri et al. 2014). Note: Spike
%           counts will be different for each group since number of cells
%           and cell firing rates are not necessarily the same.
%       F2: Slope of circular-linear regression line.
%       F3: Phase offset of circular-linear regression line.
%       F4: r^2 value of circular-linear regression line.
%
% OPTIONS:
%   See code for options.
%
% MMD
% 7/2021
% Colgin Lab

%% OPTIONS

prepForStats = 0;

degBinSz = 20; %degrees per bin for theta phase
posBinSz = 0.01; %normalized distance into place field

saveOrNot = 0; %1 to save figs, 0 to now
saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\THETA_GENERAL\phasePrecession';

cols = {'Blue', 'Red'};

%% INITIALIZE

numDegBins = 720/degBinSz; %two theta cycles
numPosBins = 1/posBinSz;

degxPos = zeros(numDegBins,numPosBins,2); %initialize, by group

ppSlopes = cell(2,1); %phase precession slopes by group
ppPhaseOffset = cell(2,1);
ppPhaseRange = cell(2,1); %pahse range = phase difference between the
% spike with the highest shifted phase and the spike with the lowest
% shifted phase (shifted phase is minus phase offset)
ppR2 = cell(2,1); %r2 values from circular regression

curDir = pwd;

%% GET DATA

for g = 1:2
    fprintf('%s\n', group(g).name)
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d\n', r, length(group(g).rat))
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day))
            for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                fprintf('\t\t\tUnit %d/%d\n', u, length(group(g).rat(r).day(d).xAllBeginUnitInfo))
                pf = group(g).rat(r).day(d).xAllBeginUnitInfo(u).pf;
                if ~isempty(pf)
                    for p = 1:length(pf) %place field ind
                        uPosxPhis = cell(2,1);
                        tetNum = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID(1);
                        
                        for b = 1:4
                            radPos = group(g).rat(r).day(d).begin(b).radPos;
                            coords = group(g).rat(r).day(d).begin(b).coords;
                            spkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                            
                            cd(group(g).rat(r).day(d).begin(b).dir)
                            lfpStruct = read_in_lfp(['CSC' num2str(tetNum) '.ncs']);
                            load(['CSC' num2str(tetNum) '_narrowThetaLfp.mat'], 'filtLfp')
                            lfpStruct.narrowThetaLfp = filtLfp;
                            
                            pp = phase_precession_circtrack_pf(spkTms, radPos, coords, pf, lfpStruct);
                            
                            if ~isempty(pp(p).spkTms)
                                uPosxPhis{1} = [uPosxPhis{1} pp(p).spkPhis];
                                uPosxPhis{2} = [uPosxPhis{2} pp(p).normSpkPos];
                                
                                tmpDegxPos = zeros(numDegBins/2,numPosBins); %initialize %first 360 degrees
                                degBins = 0:degBinSz:360;
                                
                                allPhis = pp(p).spkPhis; %get all of the spike theta phases
                                
                                for db = 1:numDegBins/2
                                    degRange = [degBins(db) degBins(db+1)];
                                    pullPos = pp(p).normSpkPos(allPhis > degRange(1) & allPhis <= degRange(2));
                                    
                                    tmpPosBins = histcounts(pullPos, 0:posBinSz:1);
                                    tmpDegxPos(db,:) = tmpPosBins;
                                end %degrees bins
                            end %if there were spikes
                            
                            tmpTwoCycle = [tmpDegxPos; tmpDegxPos];
                            degxPos(:,:,g) = degxPos(:,:,g) + tmpTwoCycle;
                            
                        end %begin
                        
                        para = circ_lin_regress(uPosxPhis{2}, deg2rad(uPosxPhis{1})); %using CZ code and methods
                        calphase = 2*pi*para(1,1)*(0:1)+para(1,2);
                        over = calphase>=2*pi;
                        under = calphase<0;
                        calphase(over) = calphase(over)-2*pi;
                        calphase(under) = calphase(under)+2*pi;
                        
                        slope = calphase(2) - calphase(1);
                        phaseOff = calphase(1);
                        phaseRange = max(deg2rad(uPosxPhis{1})) - min(deg2rad(uPosxPhis{1}));
                        R2 = circ_corrcl(deg2rad(pp(p).spkPhis), pp(p).normSpkPos);
                        
                        ppSlopes{g} = [ppSlopes{g} slope];
                        ppPhaseOffset{g} = [ppPhaseOffset{g} phaseOff];
                        ppPhaseRange{g} = [ppPhaseRange{g} phaseRange];
                        ppR2{g} = [ppR2{g} R2];
                    end %place field ind
                end %there is a place field
            end %unit
        end %day
    end %rat
end %group

keyboard

%% FIG 1 - PHASE X POS HEATMAP

figtitle = 'PhasePrecession_ThetaPhasexPostion';
figure('Name', figtitle, 'Position', [529 461 940 420])

% maxSpkCnt = max(degxPos(:));
% normDegxPos = degxPos ./ maxSpkCnt;

for g = 1:2
    subplot(1,2,g)
    imagesc(degxPos(:,:,g))
    colormap(jet)
    axis xy
    
    xlim([1 numPosBins])
    ylim([1 numDegBins])
    
    xticks(0:0.2/posBinSz:numPosBins)
    xticklabels({'0', '0.2', '0.4', '0.6', '0.8', '1'})
    xlabel('Position in place field')
    
    yticks(1:100/degBinSz:numDegBins)
    yticklabels({'0', '100', '200', '300', '400', '500', '600', '700'})
    
    if g == 1
        ylabel('Theta phase (degrees)')
    end
    
    cbr = colorbar;
    ylabel(cbr, 'Spike Counts')
    
    %     if g == 2
    %         cbr = colorbar;
    %         set(cbr, 'Position', [.93 .1 .01 .8])
    %         ylabel(cbr, 'Spike counts')
    %
    %         cbrTickInd = 50/maxSpkCnt;
    %
    %         cbr.Ticks = [0:cbrTickInd:1];
    %         cbr.TickLabels = {'0', '50', '100', '150', '200', '250', '300'};
    %     end
    
    title(group(g).name)
    
end %group

if saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 2 - SLOPE

figtitle = 'PhasePrecession_Slope';
figure('Name', figtitle, 'Position', [680 558 394 420])

jitter = 0.2;

dotplot(1:2, ppSlopes, jitter, [rgb(cols{1}); rgb(cols{2})], [0 0 0; 0 0 0]);
zero_line;

avgSlope = cellfun(@nanmean, ppSlopes);
for g = 1:2
    line([g-0.25 g+0.25], [avgSlope(g) avgSlope(g)], 'Color', 'Black', 'LineWidth', 2)
end

xticklabels({group(1).name group(2).name})
ylabel('Slope (rad/pf)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 3 - PHASE OFFSET

figtitle = 'PhasePrecession_PhaseOffset';

figure('Name', figtitle, 'Position', [680 558 394 420])

jitter = 0.2;

dotplot(1:2, ppPhaseOffset, jitter, [rgb(cols{1}); rgb(cols{2})], [0 0 0; 0 0 0]);
zero_line;

avgPhaseOff = cellfun(@nanmean, ppPhaseOffset);
for g = 1:2
    line([g-0.25 g+0.25], [avgPhaseOff(g) avgPhaseOff(g)], 'Color', 'Black', 'LineWidth', 2)
end

xticklabels({group(1).name group(2).name})
ylabel('Phase offset (rad)')


if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 4 - PHASE RANGE

figtitle = 'PhasePrecession_PhaseRange';

figure('Name', figtitle, 'Position', [680 558 394 420])

jitter = 0.2;

dotplot(1:2, ppPhaseRange, jitter, [rgb(cols{1}); rgb(cols{2})], [0 0 0; 0 0 0]);
avgPhaseRange = cellfun(@nanmean, ppPhaseRange);
for g = 1:2
    line([g-0.25 g+0.25], [avgPhaseRange(g) avgPhaseRange(g)], 'Color', 'Black', 'LineWidth', 2)
end

xticklabels({group(1).name group(2).name})
ylabel('Phase range (rad)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end


%% FIG 4 - R2

figtitle = 'PhasePrecession_R2';
figure('Name', figtitle, 'Position', [680 558 394 420])

jitter = 0.2;
dotplot(1:2, ppR2, jitter, [rgb(cols{1}); rgb(cols{2})], [0 0 0; 0 0 0]);

avgR2 = cellfun(@nanmean, ppR2);
for g = 1:2
    line([g-0.25 g+0.25], [avgR2(g) avgR2(g)], 'Color', 'Black', 'LineWidth', 2)
end

xticklabels({group(1).name group(2).name})
ylabel('R^2')


if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
    cd(curDir)
end

%% STATS

if prepForStats == 1
    keyboard
    %Watson-Williams test, as in Bieri et al 2014, done in Matlab:
    
    wwInput = zeros(2,numPosBins);
    for g = 1:2
        pullDegxPos = degxPos(1:numDegBins/2,:,g);
        [~, maxInds] = max(pullDegxPos);
        
        degBinCtrs = degBinSz/2:degBinSz:360;
        radBinCtrs = deg2rad(degBinCtrs);
        
        wwInput(g,:) = radBinCtrs(maxInds);
        
    end %group
    
    circ_wwtest(wwInput(1,:), wwInput(2,:));
    
    statSlopes = zeros(length(ppSlopes{1})+length(ppSlopes{2}),2);
    statPhaseOffset = zeros(length(ppSlopes{1})+length(ppSlopes{2}),2);
    statPhaseRange = zeros(length(ppSlopes{1})+length(ppSlopes{2}),2);
    statR2 = zeros(length(ppSlopes{1})+length(ppSlopes{2}),2);
    
    i = 0;
    for g = 1:2
        for u = 1:length(ppSlopes{g})
            i = i + 1;
                statSlopes(i,:) = [g ppSlopes{g}(u)];
                statPhaseOffset(i,:) = [g ppPhaseOffset{g}(u)];
                statPhaseRange(i,:) = [g ppPhaseRnage{g}(u)];
                statR2(i,:) = [g ppR2{g}(u)];
        end %units
    end %group
    
    
end %prep for stats

keyboard
end %function