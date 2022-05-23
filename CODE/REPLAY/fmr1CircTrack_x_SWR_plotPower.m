function fmr1CircTrack_x_SWR_plotPower(group)
% function fmr1CircTrack_x_SWR_plotPower(group)
% 
% 


%% OPTIONS

saveOrNot = 1;

minCell = 20; %number of simultanesouly recoded cells - for identifying forward or reverse events

nboot = 5000;

preWhiten = 1; %whether or not to pre-whiten the LFP (using prewhitening function in Matlab)
%% INITIALIZE
r2Thresh = 0.5; %when separating out forward vs reverse events

fRange = [2 300];
gamRange = [25 95]; %gamma
width = 7; %as in Mably et al. 2017

ripRange = [150 300];
ripInds = [ripRange(1)-(fRange(1)-1) ripRange(2)-(fRange(1)-1)];

Fs = 2000; %for lfp
order = 5; %prewhiten
useOrigNorm = 0; %see get_wavelet_power code for more information
dbConv = 0; %don't convert to db (?)

evTtls = {'Forward events', 'Reverse events'};

acEvTFR = cell(1,2); %event TFR, normalized time by group
acEvTFRbyEv = cell(2,2); %group x for/rev

PxF = cell(1,2); %power
PxFbyEv = cell(2,2); %group x ev type
pkRipFreq = cell(1,2); %for storing peak ripple frequency across event types
pkRipFreqByEv = cell(2,2);

numBins = 5000; %number of time bins when normalizing over time for evs

cols = {'Blue', 'Red'};

saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\REPLAY\TFRsharpWaveRipple';
curDir = pwd;

%% GET DATA

for g = 1:2
    fprintf('%s\n', group(g).name)
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d\n', r, length(group(g).rat))
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day))
            
            dayTFR = [];
            dayTFRbyEv = cell(1,2);
            dayPxF = [];
            dayPxFByEv = cell(1,2); %for or rev
            dayPkFreq = [];
            dayPkFreqByEv = cell(1,2);
            
            tetNums = group(g).rat(r).day(d).tetNums;
            
            for tt = 1:length(tetNums)
                fprintf('\t\t\tTetrode %d\n', tetNums(tt))
                tmpTet = tetNums(tt);
                tetTFR = [];
                tetTFRbyEv = cell(1,2);
                
                tetPxF = [];
                tetPxFByEv = cell(1,2); %for or rev
                tetPkFreq = [];
                tetPkFreqByEv = cell(1,2);
                
                for s = 2:5 %after experience
                    fprintf('\t\t\t\tSleep %d\n', s)
                    if isempty(group(g).rat(r).day(d).sleep(s).coords)
                        continue
                    end % no sleep data
                    cd(group(g).rat(r).day(d).sleep(s).dir)
                    
                    lfpStruct = read_in_lfp(['CSC' num2str(tmpTet) '.ncs']);
                    if preWhiten == 1
                        cscData = prewhitening(lfpStruct.data, order); %prewhiten
                    else
                        cscData = lfpStruct.data;
                    end %whether to prewhiten
                    
                    for i = 1:length(group(g).rat(r).day(d).sleep(s).rip)
                        if isempty(group(g).rat(r).day(d).sleep(s).rip(i).adjTms)
                            continue
                        end %no adjusted times (no spikes in detected event)
                        
                        startTm = group(g).rat(r).day(d).sleep(s).rip(i).adjTms(1);
                        endTm = group(g).rat(r).day(d).sleep(s).rip(i).adjTms(2);
                        slope = group(g).rat(r).day(d).sleep(s).rip(i).slope;
                        r2 = group(g).rat(r).day(d).sleep(s).rip(i).r2;
                        
                        adjStartInd = find(lfpStruct.ts < startTm, 1, 'Last');
                        adjEndInd = find(lfpStruct.ts < endTm, 1, 'Last');
                        
                        evLFP = cscData(adjStartInd:adjEndInd);
                        evTFR = get_wavelet_power(evLFP, Fs, fRange, width, useOrigNorm, dbConv); %get freq x time power for whole window
                        
                        tmpNormEvTFR = zeros(fRange(end)-fRange(1)+1, numBins);
                        
                        tmpTms = lfpStruct.ts(adjStartInd:adjEndInd)';
                        tmpTms = tmpTms - tmpTms(1);
                        tmpBinTm = rescale(tmpTms, 1, numBins);
                        if length(tmpBinTm) > numBins
                            keyboard
                        end
                        
                        newBinStart = 1;
                        for bn = 2:length(tmpBinTm)
                            newBinEnd = round(tmpBinTm(bn));
                            try
                                tmpNormEvTFR(:,newBinStart:newBinEnd) = repmat(evTFR(:,bn-1), [1, newBinEnd-newBinStart+1]); %tile to fit norm time distribution
                            catch; keyboard; end
                            newBinStart = newBinEnd;
                        end %bin
                        
                        tetTFR = cat(3, tetTFR, tmpNormEvTFR);
                        
                        tmpPxF = mean(evTFR, 2)';
                        tetPxF = [tetPxF; tmpPxF];
                        
                        [~, pkInd] = max(tmpPxF(ripInds(1):ripInds(2)));
                        tmpPkFreq = pkInd + ripRange(1) - 1;
                        tetPkFreq = [tetPkFreq tmpPkFreq];
                        if ~isnan(r2) && r2 > r2Thresh
                            if slope > 0 %forward event
                                sInd = 1; %positive slope
                            else %negative event
                                sInd = 2;
                            end %which direction is slope
                            tetTFRbyEv{sInd} = cat(3, tetTFRbyEv{sInd}, tmpNormEvTFR);
                            
                            tetPxFByEv{sInd} = [tetPxFByEv{sInd}; tmpPxF];
                            tetPkFreqByEv{sInd} = [tetPkFreqByEv{sInd} tmpPkFreq];
                        end %can be called a replay event
                    end %ripples - i
                end %sleep
                
                meanTetPxF = mean(tetPxF,1);
                
                tmpMean = mean(tetPxF(:));
                tmpStd = std(tetPxF(:));
                
                zMeanTetTFR = (mean(tetTFR,3) - tmpMean) / tmpStd;
                dayTFR = cat(3, dayTFR, zMeanTetTFR);
                
                tmpFunc = @(x)((mean(x,3) - tmpMean) / tmpStd); %get mean across events and zscore it
                meanTetTFRbyEv = cellfun(tmpFunc, tetTFRbyEv, 'UniformOutput', false);
                
                zTetPxF = (meanTetPxF - tmpMean) / tmpStd;
                dayPxF = [dayPxF; zTetPxF];
                dayPkFreq = [dayPkFreq mean(tetPkFreq)];
                
                tmpFunc =  @(x)((mean(x) - tmpMean) / tmpStd); %get mean and zscore it
                zTetPxFByEv = cellfun(tmpFunc, tetPxFByEv, 'UniformOutput', false);
                for sInd = 1:2
                    dayTFRbyEv{sInd} = cat(3, dayTFRbyEv{sInd}, meanTetTFRbyEv{sInd});
                    dayPxFByEv{sInd} = [dayPxFByEv{sInd}; zTetPxFByEv{sInd}];
                    dayPkFreqByEv{sInd} = [dayPkFreqByEv{sInd} mean(tetPkFreqByEv{sInd})];
                end %sInd
            end %tet
            if isempty(dayPxF)
                continue
            end %continue onto the next day
            
            acEvTFR{g} = cat(3, acEvTFR{g}, mean(dayTFR,3));
            PxF{g} = [PxF{g}; mean(dayPxF,1)];
            pkRipFreq{g} = [pkRipFreq{g} mean(dayPkFreq)];
            
            for sInd = 1:2
                acEvTFRbyEv{g,sInd} = cat(3, acEvTFRbyEv{g,sInd}, mean(dayTFRbyEv{sInd},3));
                PxFbyEv{g,sInd} = [PxFbyEv{g,sInd}; mean(dayPxFByEv{sInd},1)];
                pkRipFreqByEv{g,sInd} = [pkRipFreqByEv{g,sInd} mean(dayPkFreqByEv{sInd})];
            end %slope ind
            
        end %day
    end %rat
end %group

cd(saveDir)

keyboard
%% FIG 1 - IN EV TFR BY GROUP

tmpFunc = @(x)(mean(x,3));
meanByGroup = cellfun(tmpFunc, acEvTFR, 'UniformOutput', false);

figtitle = 'SWR_TFR_inEvent';
if preWhiten == 1
    figtitle = [figtitle '_preWhiten'];
end %add prewhiten to title
figure('Name', figtitle, 'Position', [241 558 999 420])

tmpMin = 0;
tmpMax = 0;
for g = 1:2
    subplot(1,2,g)
    
    imagesc(0:1/numBins:1, fRange(1):fRange(2), meanByGroup{g})
    axis xy
    colormap(jet)
    
    if max(meanByGroup{g}(:)) > tmpMax
        tmpMax = max(meanByGroup{g}(:));
    end
    if min(meanByGroup{g}(:)) < tmpMin
        tmpMin = min(meanByGroup{g}(:));
    end
    
    ylabel('Frequency (Hz)')
    ylim(fRange)
    xlabel('Normalized time in SWR')
    xticks([0 1])
    title(group(g).name)
    
end %group

for g = 1:2
    subplot(1,2,g)
    caxis([tmpMin tmpMax])
end

cbr = colorbar;
cbr.Position = [0.92 0.1095 0.0261 0.8167];
ylabel(cbr, 'z-scored power')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 2 - IN EV TFR BY GROUP AND EVNT DIRECTION

meanTFRbyEvType = cellfun(tmpFunc, acEvTFRbyEv, 'UniformOutput', false);

figtitle = 'SWR_TFR_byEvType';
if preWhiten == 1
    figtitle = [figtitle '_preWhiten'];
end %add prewhiten to title
figure('Name', figtitle, 'Position', [241 159 1075 819])

spMap = [1 2; 3 4];

tmpMin = inf;
tmpMax = -inf;
for g = 1:2
    for sInd = 1:2
        subplot(2,2,spMap(g,sInd))
        imagesc(0:1/numBins:1, fRange(1):fRange(2), meanTFRbyEvType{g,sInd})
        axis xy
        colormap(jet)
        
        if max(meanTFRbyEvType{g,sInd}(:)) > tmpMax
            tmpMax = max(meanTFRbyEvType{g,sInd}(:));
        end
        if min(meanTFRbyEvType{g,sInd}(:)) < tmpMin
            tmpMin = min(meanTFRbyEvType{g,sInd}(:));
        end
        
        ylabel('Frequency (Hz)')
        ylim(fRange)
        xlabel('Normalized time in SWR')
        xticks([0 1])
        title([group(g).name ' - ' evTtls{sInd}])
        
    end %slope ind
end %group

for sp = 1:4
    subplot(2,2,sp)
    caxis([tmpMin tmpMax])
end

cbr = colorbar;
cbr.Position = [0.92 0.1095 0.0261 0.8167];
ylabel(cbr, 'z-scored power')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 3 - PxF BY GROUP - ALL

figtitle = 'SWR_powerXfreq';
if preWhiten == 1
    figtitle = [figtitle '_preWhiten'];
end %add prewhiten to title
figure('Name', figtitle)
if preWhiten == 1
    figtitle = [figtitle '_preWhiten'];
end %add prewhiten to title

lh = nan(1,2);
for g = 1:2
    
    pullData = PxF{g};
    meanData = mean(pullData,1);
    tmpFunc = @(x)(mean(x,1));
    CI = bootci(nboot, {tmpFunc, pullData}, 'type', 'per');
    
    lh(g) =  plot_filled_ci(fRange(1):fRange(2), meanData, CI, rgb(cols{g}));
    
    
end %group

xlim([fRange(1) fRange(2)]);
xlabel('Frequency (Hz)')
ylabel('z-scored power')
legend(lh, {group(1).name, group(2).name}, 'Location', 'northeast')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 4 - PxF BY GROUP AND EV DIRECTION

figtitle = 'SWR_powerXfreq_byEvType';
if preWhiten == 1
    figtitle = [figtitle '_preWhiten'];
end %add prewhiten to title
figure('Name', figtitle, 'Position', [504 446 1111 471])
if preWhiten == 1
    figtitle = [figtitle '_preWhiten'];
end %add prewhiten to title

for sInd = 1:2
    subplot(1,2,sInd)
    lh = nan(1,2);
    for g = 1:2
        
        pullData = PxFbyEv{g,sInd};
        meanData = mean(pullData,1);
        tmpFunc = @(x)(mean(x,1));
        CI = bootci(nboot, {tmpFunc, pullData}, 'type', 'per');
        
        lh(g) =  plot_filled_ci(fRange(1):fRange(2), meanData, CI, rgb(cols{g}));
        
    end %group
    xlim([fRange(1) fRange(2)]);
    xlabel('Frequency (Hz)')
    ylabel('z-scored power')
    legend(lh, {group(1).name, group(2).name}, 'Location', 'northeast')
    title(evTtls{sInd})
    
end %slope

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end

%% FIG 5 - PEAK RIP FREQ

figtitle = 'SWR_peakRippleFreq';
if preWhiten == 1
    figtitle = [figtitle '_preWhiten'];
end %add prewhiten to title
figure('Name', figtitle, 'Position', [680 558 419 420])

meanPkFreq = cellfun(@mean, pkRipFreq);
semPkFreq = cellfun(@semfunct, pkRipFreq);

bgraph = bar(meanPkFreq, 'FaceColor', 'Flat');
hold on;
er = errorbar(meanPkFreq, semPkFreq);
er.Color = rgb('Black');
er.LineStyle = 'None';
for g = 1:2
    bgraph.CData(g,:) = rgb(cols{g});
end %group

xticklabels({group(1).name, group(2).name})
ylabel('Peak ripple frequency (Hz)')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end


%% FIG 6 - PEAK RIP FREQ BY EV DIR

figtitle = 'SWR_peakRippleFreq_byEvType';
if preWhiten == 1
    figtitle = [figtitle '_preWhiten'];
end %add prewhiten to title
figure('Name', figtitle)

meanPkFreq = cellfun(@mean, pkRipFreqByEv)';
semPkFreq = cellfun(@semfunct, pkRipFreqByEv)';

bgraph = bar(meanPkFreq, 'FaceColor', 'Flat');
hold on;

ngroups = 2;
nbars = 2;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for g = 1:2
    bgraph(g).CData = rgb(cols{g});
    
    x = (1:ngroups) - groupwidth/2 + (2*g-1) * groupwidth / (2*nbars);
    er = errorbar(x, meanPkFreq(:,g), semPkFreq(:,g), '.');
    er.Color = [0 0 0];
end %group

xticklabels(evTtls)
ylabel('Peak ripple frequency (Hz)')
legend({group(1).name, group(2).name}, 'Location', 'northeastoutside')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc');
    saveas(gcf, figtitle, 'fig');
    saveas(gcf, figtitle, 'png');
end


cd(curDir)

end %function