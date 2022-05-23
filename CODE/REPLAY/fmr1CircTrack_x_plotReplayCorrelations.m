function fmr1CircTrack_x_plotReplayCorrelations(group)
% function fmr1CircTrack_x_plotReplayCorrelations(group)


%% OPTIONS

saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\REPLAY\replayProperties';
saveOrNot = 0; %to save figs

downSampEvents = 1; %to randomly down sample events to to same as lowest group
if downSampEvents == 1
    rsE = RandStream('mt19937ar');
end

minFr = 1; %Hz, for cell to be included

r2Thresh = 0.5; %min r2 value for event to be included

%% INITIALIZE

% Decoding parameters
sampRate = 20000; %HzsampRate = 20000; %Hz - spike sampling rate
bayesWin = 20/1000; %20 ms time window, as in Hwaun & Colgin 2019
bayesStep = 10/1000; %10 ms time step, as in Hwaun & Colgin 2019

spatBinSz = 4; %deg

weighCorr = cell(2,1); %by group
spkTrainCorr = cell(2,1);
spkTrainCorrbyEv = cell(2,2); %by group and event direction
weighCorrbyEv = cell(2,2); %by group and event direction

slopeNames = {'Forward', 'Reverse'};

curDir = pwd; %for returning after saving figs

degBinCtrs = group(2).rat(1).day(1).binCtrs; %doesn't change across days/rats
radBinCtrs = deg2rad(degBinCtrs);

cols = {'Blue', 'Red'};

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
            
            for s = 2:length(group(g).rat(r).day(d).sleep) %after experience, so just sleep 2 on
                if ~isempty(group(g).rat(r).day(d).sleep(s).unit) %if there is anything in the sleep
                    
                    popEvents = group(g).rat(r).day(d).sleep(s).popEv; %shorten variable name
                    for i = 1:length(popEvents)
                        if isnan(popEvents(i).r2) || popEvents(i).r2 < r2Thresh
                            continue %to next event
                        end %r2 doesn't meet threshold
                        
                        if popEvents(i).slope > 0 %forward event
                            sInd = 1; %positive slope
                        else %negative event
                            sInd = 2;
                        end %which direction is slope
                        
                        pxn = popEvents(i).pxn;
                        evTms = popEvents(i).tms;
                        corrTP = replayWeighCorr(pxn, spatBinSz);
                        if ~isnan(corrTP)
                            weighCorr{g} = [weighCorr{g} abs(corrTP)];
                            weighCorrbyEv{g,sInd} = [weighCorrbyEv{g,sInd} corrTP];
                        end
                        evSpkTms = cell(1,length(uIDs));
                        uPkPos = zeros(1,length(uIDs));
                        uCntr = 0;
                        for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                            uID = group(g).rat(r).day(d).x3BeginUnitInfo(u).ID;
                            if ismember(uID, uIDs, 'row')
                                uCntr = uCntr + 1;
                                tmpEvSpks = group(g).rat(r).day(d).sleep(s).unit(u).spkTms(group(g).rat(r).day(d).sleep(s).unit(u).spkTms >= evTms(1) & group(g).rat(r).day(d).sleep(s).unit(u).spkTms <= evTms(2));
                                evSpkTms{uCntr} = tmpEvSpks;
                                
                                [~, pkInd] = max(group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap);
                                uPkPos(uCntr) = degBinCtrs(pkInd);
                            end %ismember
                        end %unit
                        
                        corrSP = replaySpkTrainCorr(evSpkTms, uPkPos);
                        if ~isnan(corrSP)
                            spkTrainCorr{g} = [spkTrainCorr{g} corrSP];
                            spkTrainCorrbyEv{g,sInd} = [spkTrainCorrbyEv{g,sInd} corrSP];
                        end
                    end %i for events
                end %if there is anything in the sleep
            end %sleep
        end %day
    end %rat
end %group

keyboard
cd(saveDir)
%% FIG 1 - WEIGHTED CORRELATION

figtitle = ['WeightedCorrelation_replay'];
if downSampEvents == 1
    figtitle = [figtitle '_downSampEvents'];
end

figure('Name', figtitle)

meanCorr = cellfun(@mean, weighCorr);
semCorr = cellfun(@semfunct, weighCorr);

bgraph = bar(meanCorr, 'FaceColor', 'Flat');

hold on;

er = errorbar(1:2, meanCorr, semCorr, semCorr, 'Color', 'Black', 'LineStyle', 'None');
xLabs = cell(1,2);
for g = 1:2
    bgraph.CData(g,:) = rgb(cols{g});
    xLabs{g} = [group(g).name ' n = ' num2str(length(weighCorr{g}))];
end

xticklabels(xLabs)
ylabel('|Weighted correlation|')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

%% FIG 2 - WEIGHTED CORR BY EV

figtitle = ['WeightedCorrelation_replay_byEv'];
if downSampEvents == 1
    figtitle = [figtitle '_downSampEvents'];
end

figure('Name', figtitle, 'Position', [478 644 793 321])

meanCorr = cellfun(@nanmean, weighCorrbyEv);
semCorr = cellfun(@nansemfunct, weighCorrbyEv);

for dirInd = 1:2
    subplot(1,2,dirInd)
    
    bgraph = bar(meanCorr(:,dirInd), 'FaceColor', 'Flat');
    
    hold on;
    
    er = errorbar(1:2, meanCorr(:,dirInd), semCorr(:,dirInd), semCorr(:,dirInd), 'Color', 'Black', 'LineStyle', 'None');
    xLabs = cell(1,2);
    for g = 1:2
        bgraph.CData(g,:) = rgb(cols{g});
        xLabs{g} = [group(g).name ' n = ' num2str(length(weighCorrbyEv{g,dirInd}))];
    end
    
    xticklabels(xLabs)
    ylabel('Weighted correlation')
    title([slopeNames{dirInd} ' Sequences'])
end %dirInd

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

%% FIG 1 - WEIGHTED CORRELATION

figtitle = ['SpikeTrainCorrelation_replay'];
if downSampEvents == 1
    figtitle = [figtitle '_downSampEvents'];
end

figure('Name', figtitle)

meanCorr = cellfun(@mean, spkTrainCorr);
semCorr = cellfun(@semfunct, spkTrainCorr);

bgraph = bar(meanCorr, 'FaceColor', 'Flat');

hold on;

er = errorbar(1:2, meanCorr, semCorr, semCorr, 'Color', 'Black', 'LineStyle', 'None');
xLabs = cell(1,2);
for g = 1:2
    bgraph.CData(g,:) = rgb(cols{g});
    xLabs{g} = [group(g).name ' n = ' num2str(length(weighCorr{g}))];
end

xticklabels(xLabs)
ylabel('Spike train correlation')

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option

%% FIG 4 - SPIKE TRAIN CORRELATION

figtitle = ['SpikeTrainCorrelation_replay_byEv'];
if downSampEvents == 1
    figtitle = [figtitle '_downSampEvents'];
end

figure('Name', figtitle, 'Position', [478 644 793 321])

meanCorr = cellfun(@nanmean, spkTrainCorrbyEv);
semCorr = cellfun(@nansemfunct, spkTrainCorrbyEv);

for dirInd = 1:2
    subplot(1,2,dirInd)
    
    bgraph = bar(meanCorr(:,dirInd), 'FaceColor', 'Flat');
    
    hold on;
    
    er = errorbar(1:2, meanCorr(:,dirInd), semCorr(:,dirInd), semCorr(:,dirInd), 'Color', 'Black', 'LineStyle', 'None');
    
    xLabs = cell(1,2);
    for g = 1:2
        bgraph.CData(g,:) = rgb(cols{g});
        xLabs{g} = [group(g).name ' n = ' num2str(length(spkTrainCorrbyEv{g,dirInd}))];
    end
    
    xticklabels(xLabs)
    ylabel('Spike train correlation')
    title([slopeNames{dirInd} ' Sequences'])
    
end %dirInd

if saveOrNot == 1
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
end %save option


end %function