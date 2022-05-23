function fmr1CircTrack_x_SWR_TFR(group)


%% OPTIONS

saveOrNot = 1;

preEvTm = 0.1; %time pre-event
postEvTm = 0.3;

minCell = 20; %number of simultanesouly recoded cells

%% INITIALIZE
r2Thresh = 0.5; %when separating out forward vs reverse events

fRange = [2 300];
gamRange = [25 95]; %gamma
width = 7; %as in Mably et al. 2017

TFRbyGroup = cell(2,1); %by group only
TFRbyEvType = cell(2,2); %group x forward vs reverse events

Fs = 2000; %for lfp
order = 5; %prewhiten
useOrigNorm = 0; %see get_wavelet_power code for more information
dbConv = 0; %don't convert to db (?)
curDir = pwd;

evTtls = {'Forward events', 'Reverse events'};
tmpRanges = [2 300; 100 300; 25 95; 55 95; 25 55; 2 25];


% cols = {'Purple', 'Green'};
cols = {'Blue', 'Red'};
rangeNames = {'All', 'Ripple', 'AllGamma', 'FastGamma', 'SlowGamma', 'ThetaDelta'};

saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\REPLAY\TFRreplayEvent';

%% GET DATA

for g = 1:2 
    fprintf('%s\n', group(g).name)
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d\n', r, length(group(g).rat))
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day))
            dayTFR = [];
            dayTFRbyEv = cell(2,1);
            
            tetNums = group(g).rat(r).day(d).tetNums;
            
            for tt = 1:length(tetNums)
                fprintf('\t\t\tTetrode %d\n', tetNums(tt))
                tmpTet = tetNums(tt);
                
                tetLFP = [];
                tetLFPbyEv = cell(2,1);
                for s = 2:5 %after experience
                    fprintf('\t\t\t\tSleep %d\n', s)
                    if isempty(group(g).rat(r).day(d).sleep(s).coords)
                        continue
                    end %if no sleep data
                    cd(group(g).rat(r).day(d).sleep(s).dir)
                    
                    lfpStruct = read_in_lfp(['CSC' num2str(tmpTet) '.ncs']);
                    cscData = prewhitening(lfpStruct.data, order); %prewhiten
                    
                    for i = 1:length(group(g).rat(r).day(d).sleep(s).rip)
                        if isempty(group(g).rat(r).day(d).sleep(s).rip(i).adjTms)
                            continue
                        end %no adjusted times (no spikes in detected event)
                        startTm = group(g).rat(r).day(d).sleep(s).rip(i).adjTms(1);
                        endTm = group(g).rat(r).day(d).sleep(s).rip(i).adjTms(2);
                        
                        if i ~= 1 && (startTm - preEvTm) <= group(g).rat(r).day(d).sleep(s).rip(i-1).tms(2) %check if another event occured in the pre-event time
                            continue
                        end %if event in pretime
                        
                        adjStartInd = find(lfpStruct.ts < startTm, 1, 'Last');
                        adjEndInd = find(lfpStruct.ts < endTm, 1, 'Last');
                        
                        pxn = group(g).rat(r).day(d).sleep(s).rip(i).pxn;
                        slope = group(g).rat(r).day(d).sleep(s).rip(i).slope;
                        
                        cutRange = [(adjStartInd - preEvTm * Fs) (adjStartInd + postEvTm * Fs)];
                        
                        if cutRange(1) < 1 || cutRange(2) > length(lfpStruct.ts) %if the times aren't within the bounds of the sleep window
                            continue
                        end %full range doesn't fit
                        
                        tmpLFP = cscData(cutRange(1):cutRange(2));
                        tetLFP = cat(2, tetLFP, tmpLFP);
                        
                        if ~isnan(group(g).rat(r).day(d).sleep(s).rip(i).r2) && group(g).rat(r).day(d).sleep(s).rip(i).r2 > r2Thresh
                            if slope > 0 %forward event
                                sInd = 1; %positive slope
                            else %negative event
                                sInd = 2;
                            end %which direction is slope
                            tetLFPbyEv{sInd} = cat(2, tetLFPbyEv{sInd}, tmpLFP);
                            
                        end %can be called a replay event
                    end %ripples - i
                end %sleep
                
                tetTFR = get_wavelet_power(tetLFP, Fs, fRange, width, useOrigNorm, dbConv); %get freq x time power for whole window
                tmpMean = mean(tetTFR(:));
                tmpStd = std(tetTFR(:));
                zTetTFR = (tetTFR - tmpMean) / tmpStd;
                dayTFR = cat(3, dayTFR, zTetTFR);
                
                tmpFunc = @(x)(get_wavelet_power(x, Fs, fRange, width, useOrigNorm, dbConv));
                tetTFRByEv = cellfun(tmpFunc, tetLFPbyEv, 'UniformOutput', false);
                
                tmpFunc = @(x)((x - tmpMean) / tmpStd);
                zTetTFRByEv = cellfun(tmpFunc, tetTFRByEv, 'UniformOutput', false);
           
                for sInd = 1:2
                    dayTFRbyEv{sInd} = cat(3, dayTFRbyEv{sInd}, zTetTFRByEv{sInd});
                end %slope ind
            end %tetNums
            if isempty(dayTFR)
                continue
            end %continue onto the next day
            meanTFR = mean(dayTFR,3);
            
            TFRbyGroup{g} = cat(3, TFRbyGroup{g}, meanTFR);
            
            tmpFunc = @(x)(mean(x,3));
            meanTFRbyEv = cellfun(tmpFunc, dayTFRbyEv, 'UniformOutput', false);
            for sInd = 1:2
                TFRbyEvType{g,sInd} = cat(3, TFRbyEvType{g,sInd}, meanTFRbyEv{sInd});
            end %slopeInd
        end %day
    end %rat
end %group
keyboard

cd(saveDir)
%% FIG 1 - BY GROUP - ALL

tmpFunc = @(x)(mean(x,3));
meanByGroup = cellfun(tmpFunc, TFRbyGroup, 'UniformOutput', false);

figtitle = 'SWR_TFR';
figure('Name', figtitle, 'Position', [241 558 999 420])

tmpMin = 0;
tmpMax = 0;
for g = 1:2
    subplot(1,2,g)
    
    imagesc(-preEvTm:1/Fs:postEvTm, fRange(1):fRange(2), meanByGroup{g})
    axis xy
    colormap(jet)
    
    if max(meanByGroup{g}(:)) > tmpMax
        tmpMax = max(meanByGroup{g}(:));
    end
    if min(meanByGroup{g}(:)) < tmpMin
        tmpMin = min(meanByGroup{g}(:));
    end
    
    line([0 0], fRange, 'Color', rgb('white'), 'LineStyle', '--')
    ylabel('Frequency (Hz)')
    ylim(fRange)
    xlabel('Time (s)')
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

%% FIG 1 - BY GROUP - GAMMA

figtitle = 'SWR_TFR_gamma';
figure('Name', figtitle, 'Position', [241 558 999 420])

tmpMin = inf;
tmpMax = -inf;
for g = 1:2
    subplot(1,2,g)
    
    imagesc(-preEvTm:1/Fs:postEvTm, gamRange(1):gamRange(2), meanByGroup{g}(gamRange(1)-(fRange(1)-1):gamRange(2)-(fRange(1)-1),:))
    axis xy
      colormap(jet)
    
      if max(max(meanByGroup{g}(gamRange(1)-(fRange(1)-1):gamRange(2)-(fRange(1)-1),:))) > tmpMax
        tmpMax = max(max(meanByGroup{g}(gamRange(1)-(fRange(1)-1):gamRange(2)-(fRange(1)-1),:)));
    end
    if min(min(meanByGroup{g}(gamRange(1)-(fRange(1)-1):gamRange(2)-(fRange(1)-1),:))) < tmpMin
        tmpMin = min(min(meanByGroup{g}(gamRange(1)-(fRange(1)-1):gamRange(2)-(fRange(1)-1),:)));
    end
    
     line([0 0], gamRange, 'Color', rgb('white'), 'LineStyle', '--')
    ylabel('Frequency (Hz)')
    ylim(gamRange)
    xlabel('Time (s)')
    
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

%% FIG 3 - BY GROUP AND EV - ALL

meanTFRbyEvType = cellfun(tmpFunc, TFRbyEvType, 'UniformOutput', false);

figtitle = 'SWR_TFR_byEvType';
figure('Name', figtitle, 'Position', [241 159 1075 819])

spMap = [1 2; 3 4];

tmpMin = inf;
tmpMax = -inf;
for g = 1:2
    for sInd = 1:2
        subplot(2,2,spMap(g,sInd))
        imagesc(-preEvTm:1/Fs:postEvTm, fRange(1):fRange(2), meanTFRbyEvType{g,sInd})
        axis xy
        colormap(jet)
        
        if max(meanTFRbyEvType{g,sInd}(:)) > tmpMax
            tmpMax = max(meanTFRbyEvType{g,sInd}(:));
        end
        if min(meanTFRbyEvType{g,sInd}(:)) < tmpMin
            tmpMin = min(meanTFRbyEvType{g,sInd}(:));
        end
        
        line([0 0], fRange, 'Color', rgb('white'), 'LineStyle', '--')
        ylabel('Frequency (Hz)')
        ylim(fRange)
        xlabel('Time (s)')
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

%% FIG 3 - BY GROUP AND EV - GAMMA

figtitle = 'SWR_TFR_byEvType_gamma';
figure('Name', figtitle, 'Position', [241 159 1075 819])

tmpMin = inf;
tmpMax = -inf;
for g = 1:2
    for sInd = 1:2
        subplot(2,2,spMap(g,sInd))
        
        tmpTFRdata = meanTFRbyEvType{g,sInd}(gamRange(1)-(fRange(1)-1):gamRange(2)-(fRange(1)-1),:);
        imagesc(-preEvTm:1/Fs:postEvTm, fRange(1):fRange(2), tmpTFRdata)
        axis xy
        colormap(jet)
        
        if max(tmpTFRdata(:)) > tmpMax
            tmpMax = max(tmpTFRdata(:));
        end
        if min(tmpTFRdata(:)) < tmpMin
            tmpMin = min(tmpTFRdata(:));
        end
        
        line([0 0], fRange, 'Color', rgb('white'), 'LineStyle', '--')
        ylabel('Frequency (Hz)')
        ylim(fRange)
        xlabel('Time (s)')
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

end %function