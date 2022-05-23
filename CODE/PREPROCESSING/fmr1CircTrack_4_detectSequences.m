function group = fmr1CircTrack_4_detectSequences(group)
% function group = fmr1CircTrack_4_detectSequences(group)
%
% PURPOSE:
%   Detect sequence events on the circle track from the continuously
%   decoded posterior probability distribution. This code uses the
%   detect_sequence_events code.
% 
% NOTE:
%   This code uses two methods allows for two methods to be used when
%   analyzing sequences from this dataset.
%       Method 1: Gets the average ratemaps from the first three run
%           sessions of the day and leaves out the last. This allows the
%           last session to serve as a separate testing dataset.
%
% MMD
% 7/2021
% Colgin Lab

%% INITIALIZE
minFr = 1; %min firing rate for unit to be included

bayesWin = 40/1000; %40 ms
bayesStep = 10/1000; %10 ms

minCellDay = 30; %minimum number of cells in a day to do the decoding

sampRate = 20000; %Hz - spike sampling rate

degBinCtrs = group(2).rat(1).day(1).binCtrs; %doesn't change across days/rats
radBinCtrs = deg2rad(degBinCtrs);

bothDir = 1;

mDesc = {'Across 3 begings', 'Across all begins for leave one out method'};


%% GET DATA

for g = 1:2
    fprintf('%s\n', group(g).name)
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d\n', r, length(group(g).rat))
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day))
            
            if length(group(g).rat(r).day(d).xAllBeginUnitInfo) < minCellDay
                fprintf('\t\t\tNot enough cells to decode\n')
                group(g).rat(r).day(d).begin(b).seq = [];
            else
                for m = 1:2
                    fprintf('\t\t\tDecoding - method %d\n', m)

                    badU = [];
                    if m == 1
                        unitInfo = group(g).rat(r).day(d).x3BeginUnitInfo;
                    else
                        unitInfo = group(g).rat(r).day(d).xAllBeginUnitInfo;
                    end
                    uIDs = zeros(length(unitInfo),2);
                    
                    rateMaps = zeros(length(unitInfo), length(unitInfo(1).smRateMap));
                    
                    for u = 1:length(unitInfo)
                        if max(unitInfo(u).rateMap)>= minFr %unit is bad if max firing rate in bin does not exceed 1
                            rateMaps(u,:) = unitInfo(u).smRateMap; %Smoothed ratemap
                            uIDs(u,:) = unitInfo(u).ID;
                        else
                            badU = [badU u]; %#ok
                        end %bad fr
                    end %units
                    rateMaps(badU,:) = [];
                    rateMaps(rateMaps==0) = 0.0001; %get rid of zeros because our Bayesian decoder can't handle 'em.
                    uIDs(badU,:) = [];
                    
                    for b = 1:4
                        fprintf('\t\t\t\tBegin %d/4\n', b)
                        group(g).rat(r).day(d).begin(b).seq(m).desc = mDesc{m};
                        radPos = group(g).rat(r).day(d).begin(b).radPos;
                        coords = group(g).rat(r).day(d).begin(b).coords;
                        
                        begStart = group(g).rat(r).day(d).begin(b).coords(1,1);
                        begEnd = group(g).rat(r).day(d).begin(b).coords(end,1);
                        begDur = begEnd - begStart;
                        
                        nTimeBins = round(begDur * sampRate); %number of bins in spike raster
                        spkRstr = zeros(size(uIDs,1), nTimeBins);
                        
                        uCntr = 0; %can't just use u since we discard some units due to firing rate (or down sampling)
                        for u = 1:length(group(g).rat(r).day(d).begin(b).unit)
                            
                            uID = unitInfo(u).ID;
                            if ismember(uID, uIDs, 'row') %if the unit wasn't discarded due to low firing ratemap
                                uCntr = uCntr + 1;
                                
                                spkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                                
                                % Take spike times and fill in the raster
                                timePassed = spkTms - begStart;
                                spkInds = round(timePassed * sampRate);
                                spkInds(spkInds<0)= []; %if spike occured before we got a time stamp for the tracking
                                
                                spkRstr(uCntr, spkInds) = 1;
                                
                            end %if unit wasn't discarded
                        end %units
                        
                        pxn = BayesianDecoder(spkRstr,rateMaps,bayesWin,bayesStep,sampRate); %Ernie's decoder
                        [seqInds, seqTms, seqSlopes] = detect_sequence_events(pxn, spkRstr, radPos, coords, radBinCtrs, bothDir);
                        
                        group(g).rat(r).day(d).begin(b).seq(m).inds = seqInds;
                        group(g).rat(r).day(d).begin(b).seq(m).tms = seqTms;
                        group(g).rat(r).day(d).begin(b).seq(m).slopes = seqSlopes;
                        group(g).rat(r).day(d).begin(b).seq(m).pxn = pxn;
                        
                    end %begin
                end %method
            end %enough cells to decode
        end %day
    end %rat
end %group

end %function