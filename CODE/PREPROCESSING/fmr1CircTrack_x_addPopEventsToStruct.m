function group = fmr1CircTrack_x_addPopEventsToStruct(group)
% function group = fmr1CircTrack_6_addReplayEventsToStruct(group)
%
% PURPOSE:
%   Detect replay events via two measures:
%       Detect population events from the binned firing rate
%       Detect SWRs from the LFP
%
% INPUT:
%   group data struct
%
% OUTPUT:
%   group data struct, with added for all beings and sleeps (with data):
%       popEvents: tms, pxn, r2, slope
%       rip: inds, tms, adjTms, pxn, r2, slope
%
% MMD
% 7/2021
% Colgin Lab

%% INITIALIZE

minFr = 1; %Hz, min firing rate for unit to be included
speedOpt = 2;

sampRate = 20000; %Hz - spike sampling rate
bayesWin = 20/1000; %20 ms time window, as in Hwaun & Colgin 2019
bayesStep = 10/1000; %10 ms time step, as in Hwaun & Colgin 2019

stdCut = 5; %standard deviations above the mean for rip power cut
durThresh = 0.5; %max dur for ripples
minDur = 0.05; %50 ms min

degBinCtrs = group(2).rat(1).day(1).binCtrs; %doesn't change across days/rats
radBinCtrs = deg2rad(degBinCtrs);

%% DO THE DATA

for g = 1:2
    fprintf('%s\n', group(g).name)
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d\n', r, length(group(g).rat))
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day))
            tetNums = group(g).rat(r).day(d).tetNums;
            if length(group(g).rat(r).day(d).xAllBeginUnitInfo) > 20 %if there are enough simultaneously recorded cells
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
                
                for b = 1:4
                    fprintf('\t\t\tBegin %d\n', b)
                    fprintf('\t\t\t\tDetecting population events\n')
                    coords = group(g).rat(r).day(d).begin(b).coords;
                    allSpkTmsbyCell = cell(1, length(group(g).rat(r).day(d).xAllBeginUnitInfo));
                    
                    badU = []; %initialize
                    for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                        if max(group(g).rat(r).day(d).xAllBeginUnitInfo(u).rateMap) >= minFr
                            allSpkTmsbyCell{1,u} = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                        else
                            allSpkTmsbyCell{1,u} = [];
                        end %reach min firing rate
                    end %unit
                    
                    popEvents = find_population_events(allSpkTmsbyCell, coords, speedOpt);
                    
                    for i = 1:size(popEvents,1)
                        startTm = popEvents(i,1);
                        endTm = popEvents(i,2);
                        
                        nEvBins = round((endTm-startTm) * sampRate); %number of bins in spike raster
                        spkRstr = zeros(size(uIDs,1), nEvBins);
                        timeaxis = 0:bayesStep:(endTm-startTm);
                        
                        uCntr = 0; %can't just use u since we discard some units due to firing rate (or down sampling)
                        for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                            
                            uID = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID;
                            if ismember(uID, uIDs, 'row') %if the unit wasn't discarded due to low firing ratemap
                                uCntr = uCntr + 1;
                                allSpkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                                tmpEvSpkTms = allSpkTms(allSpkTms>=startTm & allSpkTms<=endTm);
                                
                                timePassed = tmpEvSpkTms - startTm;
                                spkInds = round(timePassed * sampRate);
                                spkInds(spkInds==0)=1;
                                spkRstr(uCntr, spkInds) = 1;
                            end %if unit wasn't discarded
                        end %units
                        
                        pxn = BayesianDecoder(spkRstr,rateMaps,bayesWin,bayesStep,sampRate); %Ernie's decoder
                        [nWin, winStartInds] = find_num_windows(size(spkRstr,2), bayesWin*sampRate, bayesStep*sampRate); %number of win can be miscalculated by decoder
                        if winStartInds(end)+bayesWin*sampRate < size(spkRstr,2)
                            nWin = nWin + 1;
                            winStartInds(end+1) = winStartInds(end)+bayesStep*sampRate; %#ok
                        end %do we need to add another to nWin
                        pxn = pxn(:,1:nWin);
                        
                        maxVals= max(pxn);
                        bins2use = 1:size(pxn,2);
                        bins2use(isnan(maxVals)) = []; %don't use these - not better than chance
                        
                        if length(bins2use)/size(pxn,2) >= .8 %as long as spikes in at least 80% of bins - I don't think find_population_events does this as of writing this code but it should
                            [r2, ~, ~, ~, slope] = Cir_reg(pxn, radBinCtrs, timeaxis, bins2use); %CZ code - slope in rad/s
                        else
                            r2 = nan;
                            slope = nan;
                        end %if at least 80% of bins have spikes
                        
                        group(g).rat(r).day(d).begin(b).popEv(i).tms = popEvents(i,:);
                        group(g).rat(r).day(d).begin(b).popEv(i).pxn = pxn;
                        group(g).rat(r).day(d).begin(b).popEv(i).r2 = r2;
                        group(g).rat(r).day(d).begin(b).popEv(i).slope = slope;
                    end %pop events
                    
                    fprintf('\t\t\t\tDetecting SWR events\n')
                    ripInds = [];
                    ripTms = [];
                    
                    cd(group(g).rat(r).day(d).begin(b).dir)
                    CA1eeg = {};
                    for t = 1:length(tetNums)
                        CSCname = ['CSC' num2str(tetNums(t)) '.ncs'];
                        lfpStruct = read_in_lfp(CSCname);
                        CA1eeg = cat(1, CA1eeg, lfpStruct.data);
                    end
                    
                    [ripOnInds, ripOffInds] = DetectRipples_v4(CA1eeg,lfpStruct.Fs, stdCut, durThresh); %fs same for all tet
                    ripInds(:,1) = ripOnInds;
                    ripInds(:,2) = ripOffInds;
                    ripTms(:,1) = lfpStruct.ts(ripOnInds);%ts same for all tets
                    ripTms(:,2) = lfpStruct.ts(ripOffInds);
                    
                    for i = 1:size(ripOnInds,1)
                        group(g).rat(r).day(d).begin(b).rip(i).inds = ripInds(i,:);
                        group(g).rat(r).day(d).begin(b).rip(i).tms = ripTms(i,:);
                        
                        startTm = ripTms(i,1);
                        endTm = ripTms(i,2);
                        nEvBins = round((endTm-startTm) * sampRate); %number of bins in spike raster
                        spkRstr = zeros(size(uIDs,1), nEvBins);
                        
                        uCntr = 0; %can't just use u since we discard some units due to firing rate (or down sampling)
                        for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                            uID = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID;
                            if ismember(uID, uIDs, 'row') %if the unit wasn't discarded due to low firing ratemap
                                uCntr = uCntr + 1;
                                
                                allSpkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                                tmpEvSpkTms = allSpkTms(allSpkTms>=startTm & allSpkTms<=endTm);
                                timePassed = tmpEvSpkTms - startTm;
                                spkInds = round(timePassed * sampRate);
                                spkInds(spkInds==0)=1;
                                
                                spkRstr(uCntr, spkInds) = 1;
                            end %if unit wasn't discarded
                        end %units
                        
                        pxn = BayesianDecoder(spkRstr,rateMaps,bayesWin,bayesStep,sampRate); %Ernie's decoder
                        [nWin, winStartInds] = find_num_windows(size(spkRstr,2), bayesWin*sampRate, bayesStep*sampRate); %number of win can be miscalculated by decoder
                        if winStartInds(end)+bayesWin*sampRate < size(spkRstr,2)
                            nWin = nWin + 1;
                            winStartInds(end+1) = winStartInds(end)+bayesStep*sampRate; %#ok
                        end
                        winStartTms = startTm + winStartInds/sampRate;
                        winEndTms = winStartTms + bayesWin;
                        pxn = pxn(:,1:nWin);
                        
                        maxVals = max(pxn);
                        goodBins = find(~isnan(maxVals));
                        if isempty(goodBins) || round((winEndTms(goodBins(end))-winStartTms(goodBins(1))), 2) < minDur
                            group(g).rat(r).day(d).begin(b).rip(i).pxn = [];
                            group(g).rat(r).day(d).begin(b).rip(i).r2 = nan;
                            continue
                        else
                            pxn = pxn(:,goodBins(1):goodBins(end)); %cut pxn so first and last bins have spikes
                            
                            maxVals = max(pxn);
                            
                            group(g).rat(r).day(d).begin(b).rip(i).adjTms = [winStartTms(goodBins(1)) winEndTms(goodBins(end))]; %get adjusted time bins
                            
                            bins2use = 1:size(pxn,2);
                            bins2use(isnan(maxVals)) = []; %don't use these - no spikes in this bin
                            timeaxis = 0:bayesStep:(endTm-startTm); %for calculating slope
                            
                            if length(bins2use)/size(pxn,2) >= .8 %as long as spikes in at least 80% of bins
                                [r2, ~, ~, ~, slope] = Cir_reg(pxn, radBinCtrs, timeaxis, bins2use); %CZ code - slope in rad/s
                            else
                                r2 = nan;
                                slope = nan;
                            end %80% of bins have spikes
                        end %there's anything in the pxn
                        group(g).rat(r).day(d).begin(b).rip(i).pxn = pxn;
                        group(g).rat(r).day(d).begin(b).rip(i).r2 = r2;
                        group(g).rat(r).day(d).begin(b).rip(i).slope = slope;
                    end %swr events
                end %begin
                
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~SLEEPS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                for s = 1:5
                    
                    if ~isempty(group(g).rat(r).day(d).sleep(s).unit)
                        fprintf('\t\t\tSleep %d\n', s)
                        fprintf('\t\t\t\tDetecting population events\n')
                        coords = group(g).rat(r).day(d).sleep(s).coords;
                        allSpkTmsbyCell = cell(1, length(group(g).rat(r).day(d).xAllBeginUnitInfo));
                        
                        badU = []; %initialize
                        for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                            if max(group(g).rat(r).day(d).xAllBeginUnitInfo(u).rateMap) >= minFr
                                allSpkTmsbyCell{1,u} = group(g).rat(r).day(d).sleep(s).unit(u).spkTms;
                            else
                                allSpkTmsbyCell{1,u} = [];
                            end %reach min firing rate
                        end %unit
                        
                        popEvents = find_population_events(allSpkTmsbyCell, coords, speedOpt);
                        
                        for i = 1:length(popEvents)
                            startTm = popEvents(i,1);
                            endTm = popEvents(i,2);
                         
                            nEvBins = round((endTm-startTm) * sampRate); %number of bins in spike raster
                            spkRstr = zeros(size(uIDs,1), nEvBins);
                            timeaxis = 0:bayesStep:(endTm-startTm);
                            
                            uCntr = 0; %can't just use u since we discard some units due to firing rate (or down sampling)
                            for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                                uID = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID;
                                if ismember(uID, uIDs, 'row') %if the unit wasn't discarded due to low firing ratemap
                                    uCntr = uCntr + 1;
                                    allSpkTms = group(g).rat(r).day(d).sleep(s).unit(u).spkTms;
                                    tmpEvSpkTms = allSpkTms(allSpkTms>=startTm & allSpkTms<=endTm);
                                    
                                    timePassed = tmpEvSpkTms - startTm;
                                    spkInds = round(timePassed * sampRate);
                                    spkInds(spkInds==0)=1;
                                    spkRstr(uCntr, spkInds) = 1;
                                end %if unit wasn't discarded
                            end %units
                            pxn = BayesianDecoder(spkRstr,rateMaps,bayesWin,bayesStep,sampRate); %Ernie's decoder
                            [nWin, winStartInds] = find_num_windows(size(spkRstr,2), bayesWin*sampRate, bayesStep*sampRate); %number of win can be miscalculated by decoder
                            if winStartInds(end)+bayesWin*sampRate < size(spkRstr,2)
                                nWin = nWin + 1;
                                winStartInds(end+1) = winStartInds(end)+bayesStep*sampRate; %#ok
                            end %do we need to add another to nWin
                            pxn = pxn(:,1:nWin);
                            
                            maxVals= max(pxn);
                            bins2use = 1:size(pxn,2);
                            bins2use(isnan(maxVals)) = []; %don't use these - not better than chance
                            
                            if length(bins2use)/size(pxn,2) >= .8 %as long as spikes in at least 80% of bins - I don't think find_population_events does this as of writing this code but it should
                                 [r2, ~, ~, ~, slope]   = Cir_reg(pxn, radBinCtrs, timeaxis, bins2use); %CZ code - slope in rad/s
                            else
                                r2 = nan;
                            end %if at least 80% of bins have spikes
                            
                            group(g).rat(r).day(d).sleep(s).popEv(i).tms = popEvents(i,:);
                            group(g).rat(r).day(d).sleep(s).popEv(i).pxn = pxn;
                            group(g).rat(r).day(d).sleep(s).popEv(i).r2 = r2;
                            group(g).rat(r).day(d).sleep(s).popEv(i).slope = slope;
                        end %pop events
                        
                        fprintf('\t\t\t\tDetecting SWR events\n')
                        ripInds = [];
                        ripTms = [];
                        
                        cd(group(g).rat(r).day(d).sleep(s).dir)
                        CA1eeg = {};
                        for t = 1:length(tetNums)
                            CSCname = ['CSC' num2str(tetNums(t)) '.ncs'];
                            lfpStruct = read_in_lfp(CSCname);
                            CA1eeg = cat(1, CA1eeg, lfpStruct.data);
                        end
                        
                        [ripOnInds, ripOffInds] = DetectRipples_v4(CA1eeg,lfpStruct.Fs, stdCut, durThresh); %fs same for all tet
                        ripInds(:,1) = ripOnInds;
                        ripInds(:,2) = ripOffInds;
                        ripTms(:,1) = lfpStruct.ts(ripOnInds);%ts same for all tets
                        ripTms(:,2) = lfpStruct.ts(ripOffInds);
                        
                        for i = 1:length(ripOnInds)
                            group(g).rat(r).day(d).sleep(s).rip(i).inds = ripInds(i,:);
                            group(g).rat(r).day(d).sleep(s).rip(i).tms = ripTms(i,:);
                            
                            startTm = ripTms(i,1);
                            endTm = ripTms(i,2);
                            nEvBins = round((endTm-startTm) * sampRate); %number of bins in spike raster
                            spkRstr = zeros(size(uIDs,1), nEvBins);
                            
                            uCntr = 0; %can't just use u since we discard some units due to firing rate (or down sampling)
                            for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                                uID = group(g).rat(r).day(d).xAllBeginUnitInfo(u).ID;
                                if ismember(uID, uIDs, 'row') %if the unit wasn't discarded due to low firing ratemap
                                    uCntr = uCntr + 1;
                                    
                                    allSpkTms = group(g).rat(r).day(d).sleep(s).unit(u).spkTms;
                                    tmpEvSpkTms = allSpkTms(allSpkTms>=startTm & allSpkTms<=endTm);
                                    timePassed = tmpEvSpkTms - startTm;
                                    spkInds = round(timePassed * sampRate);
                                    spkInds(spkInds==0)=1;
                                    
                                    spkRstr(uCntr, spkInds) = 1;
                                end %if unit wasn't discarded
                            end %units
                            
                            pxn = BayesianDecoder(spkRstr,rateMaps,bayesWin,bayesStep,sampRate); %Ernie's decoder
                            [nWin, winStartInds] = find_num_windows(size(spkRstr,2), bayesWin*sampRate, bayesStep*sampRate); %number of win can be miscalculated by decoder
                            if winStartInds(end)+bayesWin*sampRate < size(spkRstr,2)
                                nWin = nWin + 1;
                                winStartInds(end+1) = winStartInds(end)+bayesStep*sampRate; %#ok
                            end
                            winStartTms = startTm + winStartInds/sampRate;
                            winEndTms = winStartTms + bayesWin;
                            pxn = pxn(:,1:nWin);
                            
                            maxVals = max(pxn);
                            goodBins = find(~isnan(maxVals));
                            if isempty(goodBins) || round((winEndTms(goodBins(end))-winStartTms(goodBins(1))), 2) < minDur
                                group(g).rat(r).day(d).sleep(s).rip(i).pxn = [];
                                group(g).rat(r).day(d).sleep(s).rip(i).r2 = nan;
                                continue
                            else
                                pxn = pxn(:,goodBins(1):goodBins(end)); %cut pxn so first and last bins have spikes
                                maxVals = max(pxn);
                                
                                group(g).rat(r).day(d).sleep(s).rip(i).adjTms = [winStartTms(goodBins(1)) winEndTms(goodBins(end))]; %get adjusted time bins
                                
                                bins2use = 1:size(pxn,2);
                                bins2use(isnan(maxVals)) = []; %don't use these - no spikes in this bin
                                timeaxis = 0:bayesStep:(endTm-startTm); %for calculating slope
                                
                                if length(bins2use)/size(pxn,2) >= .8 %as long as spikes in at least 80% of bins
                                    [r2, ~, ~, ~, slope] = Cir_reg(pxn, radBinCtrs, timeaxis, bins2use); %CZ code - slope in rad/s
                                else
                                    r2 = nan;
                                    slope = nan;
                                end %80% of bins have spikes
                                
                                group(g).rat(r).day(d).sleep(s).rip(i).pxn = pxn;
                                group(g).rat(r).day(d).sleep(s).rip(i).r2 = r2;
                                group(g).rat(r).day(d).sleep(s).rip(i).slope = slope;
                            end %there were any good bins
                        end %swr events
                    end %if there are sleep cells to use
                end %sleeps
            end %if there are enough cells
        end %day
    end %rat
end %group

end %function