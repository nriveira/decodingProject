function [seqInds, seqTms, seqSlopes] = detect_sequence_events(pxn, spkRstr, radPos, coords, radBinCtrs, bothDir)
% function [seqInds, seqTms, seqSlopes] = detect_sequence_events(pxn, spkRstr, radPos, coords, radBinCtrs, bothDir)
%
% PURPOSE:
%   To detect sequence events using the method stated in Zheng et al. 2021.
%       Works for circle track data, including free runs where the rat
%       crosses 0 degrees (unlike in DMP task). For code reference, see
%       DetectSequenceEvents_cz. NOTE: This code institutes an additional
%       check that the rat was traveling at least 5 cm/s during the
%       sequence time.
%
% INPUTS:
%   pxn = output from BayesianDecoder function
%   spkRstr = SAME spike raster used as input for Bayesian Decoder function
%       to get pxn
%   radPos = radial position of the rat over the same time period as pxn
%   coords = coordinates for the rat over same time period as pxn
%   radBinCtrs = radial position of each bin for pxn
%   bothDir = 1 for both forward and reverse sequences, 0 for forward only
%
% OUTPUTS:
%	seqTimes: (:,1) = onset time, (:,2) = offset time for sequence
%       if bothDir == 1, cell with {1} = forward sequences, {2} = reverse
%           sequences
%       if bothDir == 0, array with times for forward sequences only
%   seqSlopes: slope in rad/s for each event
%       if bothDir == 1, cell with {1} = slopes for forward sequences, {2}
%           = slopes for reverse sequences
%       if bothDir == 0, array with slopes for forward sequences only
%
% OPTIONS:
%   CHECK FUNCTION FOR THRESHOLDS! Set to match those set in Zheng et al.
%       2021.
%
% MMD
% 7/2021
% Colgin Lab

%% OPTIONS/INITIALIZE

maxJumpThr = 1.4; %rad - "estimated positions between adjacent time bins did not exceed 1.4 rad"
timeWin = 6; %6 subsequent time bins with spike 
timeStep = 1; %bin step
distanceThr = 0.07; %"distance between first and last estimated position within seqeunce was more or equal to 0.07 rad"
% jumpMeanThr = 0.5;
runThresh = 5; %cm/s - THIS WAS NOT ORIGINALLY USED IN ZHENG ET AL. 2021

bayesWin = 40/1000; %40 ms
bayesStep = 10/1000; %10 ms
sampRate = 20000; %Hz - spike sampling rate for spike raster input to Bayesian decoder

minSpk = 5; %5 spikes in an event
minCell = 3; %3 cells must spike in event

ppToFitThr = 0.35; %rad, "at least 60% of  the total post prob needed to be no mroe than 0.35 rad away from the fitted traj line"
propFitThr = 0.6;

fitToActThr = 0.35; %rad, minimal distance between the fitted trajectory and rats' true postion had to be less than or equal to 0.35 rad

if max(radPos(:,2)) > 2*pi
    radPos(:,2) = deg2rad(radPos(:,2));
end

if bothDir == 1
    seqTimes = cell(2,1);
    seqSlopes = cell(2,1);
    dirs = [1 2];
else
    seqTimes = [];
    seqSlopes = [];
    dirs = 1;
end %both dir or not

onset = cell(2,1);
offset = cell(2,1);

%% IDENTIFY POTENTIAL EVENTS FROM PXN

for ii = 1:timeStep:size(pxn,2)
    range = ii:ii+timeWin-1;
    if max(range)>size(pxn,2)
        break
    end %min range fits in pxn
    
    emptyBin = isnan(pxn(1,range));
    if any(emptyBin)
        continue %to next ind
    end %any nan
    
    [~, decodedPos] = max(pxn(:,range)); %decoded position based on max likely
    radDecPos = radBinCtrs(decodedPos);

    jump = nan(1,length(radDecPos)-1);
    for ij = 1:length(jump)
        jump(ij) = circ_dist(radDecPos(ij+1), radDecPos(ij));
    end %get jump
    
    maxJump = nanmax(abs(jump));
    
    seqDist = circ_dist(radDecPos(end), radDecPos(1));
    distFromStart = circ_dist(radDecPos, radDecPos(1));
    
    if bothDir == 1
        if maxJump <= maxJumpThr && abs(seqDist) >= distanceThr  
            if seqDist > 0
                dirInd = 1;
            else
                dirInd = 2;
            end %direction of seq - potential forward or rev?
            
            if dirInd == 1 && min(distFromStart) >= 0 || dirInd == 2 && max(distFromStart) <=0
                onset{dirInd} = cat(1,onset{dirInd},ii);
                offset{dirInd} = cat(1,offset{dirInd},range(end));
            end
        end %qualifies
    else
        dirInd = 1;
        if maxJump <= maxJumpThr && seqDist >= distanceThr && min(distFromStart) >= 0
            %                 seq = [seq; onOff];
            onset{dirInd} = cat(1,onset{dirInd},ii);
            offset{dirInd} = cat(1,offset{dirInd},range(end));
        end %meets crit
    end %both directions
end %pxn bins

%% COMBINE OVERLAPPING EVENTS

for dirInd = dirs %not combining forward events with reverse events
    if length(onset{dirInd}) > 1
        isi = onset{dirInd}(2:end)-offset{dirInd}(1:end-1);
        merge = find(isi <= 0);
        onset{dirInd}(merge+1) = [];
        offset{dirInd}(merge) = [];
    end
end %dir ind

%% CHECK CELLS/SPIKES

%make nsameple a odd number - done in BayesianDecoder function
if mod(sampRate,2) == 0
    nsample = bayesWin*sampRate+1;
else
    nsample = bayesWin*sampRate;
end

tRange = floor(nsample/2)+1:bayesStep*sampRate:size(spkRstr,2)-floor(nsample/2); %time range
for dirInd = dirs
    badInds = []; %initialize
    
    for s = 1:length(onset{dirInd})
        rangeOn = tRange(onset{dirInd}(s))-floor(nsample/2);
        try
            rangeOff = tRange(offset{dirInd}(s))+floor(nsample/2)-1;
        catch
            if offset{dirInd}(s) > length(tRange)
                rangeOff = size(spkRstr,2);
            end
        end %time range doesn't fit
        
        pullRstr = spkRstr(:,rangeOn:rangeOff);
        
        numSpksByCell = sum(pullRstr,2);
        
        if sum(pullRstr(:)) < minSpk || length(numSpksByCell(numSpksByCell~=0)) < minCell %if not good
            badInds = [badInds; s];
        end %meet spike and cell crit
    end %all seq
    
    if ~isempty(onset{dirInd})
        onset{dirInd}(badInds) = [];
        offset{dirInd}(badInds) = [];
    end %anything to delete
end %dirInd

%% GET TIMES

runSpd = get_runspeed(coords);
smRunSpd = smooth_runspeed(runSpd);

onsetTm = cell(2,1);
offsetTm = cell(2,1);

[nWin, winStartInds] = find_num_windows(size(spkRstr,2), bayesWin*sampRate, bayesStep*sampRate);
if winStartInds(end)+bayesWin*sampRate < size(spkRstr,2)
    nWin = nWin + 1;
    winStartInds(end+1) = winStartInds(end)+bayesStep*sampRate; %okay
end
winStartTms = radPos(1,1) + winStartInds/sampRate;
winEndTms = winStartTms + bayesWin;

for dirInd = dirs
    badInds = []; %initialize
    
    for s = 1:length(onset{dirInd})
        if offset{dirInd}(s) > nWin
            badInds = [badInds s];
        else
            tmpOn = onset{dirInd}(s);
            tmpOff = offset{dirInd}(s);
            
            onTm = winStartTms(tmpOn);
            offTm = winEndTms(tmpOff);
            onsetTm{dirInd} = [onsetTm{dirInd} onTm];
            offsetTm{dirInd} = [offsetTm{dirInd} offTm];         
        end %within actual nWin
    end %seq
    
    onset{dirInd}(badInds) = [];
    offset{dirInd}(badInds) = [];
end %dir ind

%% CHEK RUN SPEED

for dirInd = dirs
    badInds = []; %initialize
    for s = 1:length(onset{dirInd})
         onTm = onsetTm{dirInd}(s);
         offTm = offsetTm{dirInd}(s);
        
        startInd = find(smRunSpd(:,1) <= onTm, 1, 'Last');
        endInd = find(smRunSpd(:,1) <= offTm, 1, 'Last');
        
        seqSpd = smRunSpd(startInd:endInd, 2);
        if min(seqSpd) < runThresh %if speed drops below threshold during the sequence
            badInds = [badInds s];
        end %speed check
    end %sequences

    onset{dirInd}(badInds) = [];
    onsetTm{dirInd}(badInds) = [];
    offset{dirInd}(badInds) = [];
    offsetTm{dirInd}(badInds) = [];
    
end %dirInd

%% CHECK SUFFICIENT % PXN IS CLOSE TO FITTED LINE / CHECK FITTED CLOSE ENOUGH TO TRUE POS

seqSlopes = cell(2,1); %initialize
for dirInd = dirs
    badInds = []; %re-initialize

    for s = 1:length(onset{dirInd})
        
        tmpOn = onset{dirInd}(s);
        tmpOff = offset{dirInd}(s);
        
        pullPxn = pxn(:,tmpOn:tmpOff);
        
        [maxVals, decodedPos] = max(pullPxn);
        timeVals = 0:bayesStep:(size(pullPxn,2)*bayesStep)-bayesStep;
        
        decodedPos(maxVals==1/size(spkRstr,1)) = []; %no better than chance
        timeVals(maxVals==1/size(spkRstr,1)) = [];
        radDecPos = radBinCtrs(decodedPos);
        
        para = circ_lin_regress(timeVals, radDecPos, 8); %bound is 8 for determining slope in cz code
        
        calphase = 2*pi*para(1,1)*timeVals + para(1,2);
        slope = para(1,1)*2*pi;
        seqSlopes{dirInd} = [seqSlopes{dirInd}; slope];
        
        if dirInd == 1 && slope < 0
            badInds = [badInds s];
            continue
        elseif dirInd == 2 && slope > 0
            badInds = [badInds s];
            continue
        end %slope facing right way
        
        over = calphase>=2*pi;
        under = calphase<0;
        calphase(over) = calphase(over)-2*pi;
        calphase(under) = calphase(under)+2*pi;
        
%         fitDist = zeros(1,length(timeVals));
        sumPxn = 0;
        for x = 1:length(timeVals)
            fitPos = calphase(x); %fit position
            
            fitPosMin = fitPos - ppToFitThr;
            fitPosMax = fitPos + ppToFitThr;
            
            [~,minInd] = min(abs(circ_dist(radBinCtrs, fitPosMin)-0));
            [~,maxInd] = min(abs(circ_dist(radBinCtrs, fitPosMax)-0));
            
            if maxInd > minInd
                sumPxn = sumPxn + sum(pullPxn(minInd:maxInd,x));
            else
                sumPxn = sumPxn + sum(pullPxn(minInd:end,x)) + sum(pullPxn(1:maxInd,x)); %fit crosses 0
            end 
            
        end %time vals
        
        sumAll = sum(pullPxn(:));
        propWithDist = sumPxn/sumAll; %proportion of pxn within distance
        
        if propWithDist < propFitThr
            badInds = [badInds s];
            continue
        end %not high enough %
        
        %     keyboard
        %     figure;
        %     plot(timeVals, radDecPos)
        % %     y_reg = 2*pi*para(1)*timeVals + para(2);
        %     hold on;
        %     plot(timeVals, calphase)
        %     ylim([0 2*pi])
     
         onTm = onsetTm{dirInd}(s); %time seq starts in s
         offTm = offsetTm{dirInd}(s);
        
        startInd = find(radPos(:,1) <= onTm, 1, 'Last');
        endInd = find(radPos(:,1) <= offTm, 1, 'Last');
        
        actPos = mean(radPos(startInd:endInd,2));
        
        minDist = min(abs(circ_dist(calphase, actPos)));
        
        if minDist > fitToActThr
            badInds = [badInds s];
        end %doesn't meet min dist crit
    end %seq
    
    onset{dirInd}(badInds) = [];
    onsetTm{dirInd}(badInds) = [];
    offset{dirInd}(badInds) = [];
    offsetTm{dirInd}(badInds) = [];
    seqSlopes{dirInd}(badInds) = [];
    
end %dirInd

%% PREPARE OUTPUT

if bothDir == 1
    for dirInd = dirs
        seqInds{dirInd} = zeros(length(onset{dirInd}),2);
        seqInds{dirInd}(:,1) = onset{dirInd};
        seqInds{dirInd}(:,2) = offset{dirInd};
        
        seqTms{dirInd} = zeros(length(onset{dirInd}),2);
        seqTms{dirInd}(:,1) = onsetTm{dirInd};
        seqTms{dirInd}(:,2) = offsetTm{dirInd};
    end %direction
else
    dirInd = 1;
    
    seqInds = zeros(length(onset{dirInd}),2);
    seqInds(:,1) = onset{dirInd};
    seqInds(:,2) = offset{dirInd};
    
     seqTms = zeros(length(onset{dirInd}),2);
    seqTms(:,1) = onsetTm{dirInd};
    seqTms(:,2) = offsetTm{dirInd};
    
    seqSlopes = seqSlopes{dirInd};
    
end %both directions

end %function