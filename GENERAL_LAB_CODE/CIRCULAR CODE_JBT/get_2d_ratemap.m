function [rateMap, spkCnt, timePerBin] = get_2d_ratemap(spkTms, coords, xBnds, yBnds, spatBinSz, plotOrNot, velFilt, durCrit)
% function [rateMap, spkCnt, timePerBin] = get_2d_ratemap(spkTms, coords, xBnds, yBnds, spatBinSz, plotOrNot, velFilt, durCrit)
%
% PURPOSE:
%   To get the rate-map matrix for a cell.
%
% INPUT:
%      spkTms = spike time vector (n x 1) in seconds
%      coords = coordinate matrix (n x 3) where (:,1) = frametimes, (:,2) = xPos, and (:,3) = yPos
%               **xPos and yPos MUST be positive
%      xBnds = 1x2 vector indicating spatial boundaries along x dimension
%      yBnds = 1x2 vector indicating spatial boundaries along y dimension
%              - These two 'bnds' variables were added to allow for a ratemap to have a span around zero
%                e.g., xBnds = [-60 60]; yBnds = [-60 60]; see usage in enrichProj_1...
%   spatBinSz = size, in cm, of spatial bins.
%   plotOrNot = optional binary input indicating whether(1) or not(0) to plot the results
%                - if not provided, function defaults to 0 (do not plot)
%     velFilt = optional binary input indicating whether(1) or not(0) to only consider spikes occuring
%               when the rat was moving >= 5 cm/s
%                - if not provided, defaults to 0 (do not velocity filter)
%     durCrit = optional binary indicating whether(1) or not(0) to only consider spatial bins that the
%               rat visited for at least a specified amount of time (150 ms per Schlesiger et al., 2015)
%                - The purpose of this is not to over-inflate bins due to low occupancy
%                - If not provided, default is 0 (do not enforce duration criteria)
%
% OUTPUT:
%     rateMap = raw, unsmoothed ratemap
%      spkCnt = spike counts for each spatial bin
%  timePerBin = time, in s, per spatial bin
%
% JB Trimper
% 7/19/17
% Colgin Lab


minBinDur = 150; %ms -- minimum time spent in a bin for the firing rate there to be calculated
%                       *Employed only if durCrit == 1


%% DEFAULT, IF NECESSARY, TO NOT PLOTTING
if nargin < 5 || ~exist('plotOrNot', 'var')
    plotOrNot = 0;
end


%% DEFAULT, IF NECESSARY, TO NOT VELOCITY FILTERING
if nargin < 6 || ~exist('velFilt', 'var')
    velFilt = 0;
end

%% DEFAULT, IF NECESSARY, TO NOT COMPARING TO DURATION CRITERIA
if nargin < 7 || ~exist('durCrit', 'var')
    durCrit = 0;
end


%% GET THE SPATIAL BINS BASED ON THE BIN SIZE AND THE ARENA DIMENSIONS
numXBins = floor(diff(xBnds)/spatBinSz); % Number of spatial bins is rounded down range of x or y values...
numYBins = floor(diff(yBnds)/spatBinSz); % divided by the spatial bin size

xBinEdges = xBnds(1):spatBinSz:xBnds(2);
yBinEdges = yBnds(1):spatBinSz:yBnds(2);



%% PLOT THE RAT'S PATH
if plotOrNot == 1
    figure('Position', [411 564 1109 389])
    subplot(1,2,1);
    plot(coords(:,2), coords(:,3), 'Color', [0.5 0.5 0.5]);
    hold on;
    spkPos = zeros(2,length(spkTms));
    for st = 1:length(spkTms)
        tmpTm = spkTms(st);
        posInd = find(coords(:,1)<=tmpTm, 1, 'Last');
        spkPos(:,st) = coords(posInd,2:3);
    end
    plot(spkPos(1,:), spkPos(2,:), '.r')
    title('Path');
    xlim([xBnds(1) xBnds(2)]);
    xlim([yBnds(1) yBnds(2)]);
    xlabel('Position (cm)')
    ylabel('Position (cm)')
    set(gca, 'FontSize', 14);
    axis square;
end



%% IF DOING VELOCITY FILTERED RATEMAP, GET A BINARY THAT INDICATES WHEN THE RAT WAS MOVING
if velFilt == 1
    
    runThresh = 5; %cm/s -- threshold for movement
    
    instRs = get_runspeed(coords);
    smRs = smooth_runspeed(instRs);
    velFiltBnry = zeros(1,length(smRs));
    velFiltBnry(smRs(:,2)>=runThresh) = 1;
    
end



% This is used for making the time-stamps for when the rat was in each spatial bin more precise
halfFt = mean(diff(coords(:,1)))/2; %half of the frame-rate
%      Mark start time for when the rat was in each bin as the frame-time (mid-point) - half the frame rate
%      and mark the end time as frame-time (mid-point) + half the frame rate



%% GET THE TIMES WHEN THE RAT WAS IN EACH SPATIAL BIN
binTimes = cell(numXBins, numYBins); %catcher for time-stamps when rat is in each bin
binTimeSum = zeros(numXBins, numYBins);
for x = 1:numXBins
    xRange = [xBinEdges(x) xBinEdges(x+1)];
    for y = 1:numYBins
        yRange = [yBinEdges(y) yBinEdges(y+1)];
        
        binInds = coords(:,2)>=xRange(1) & coords(:,2)<xRange(2) & coords(:,3)>=yRange(1) & coords(:,3)<yRange(2);
        %                   Indices for when the rat is in each spatial bin
        
        tmpTimeBinary = zeros(1,length(coords(:,1))); %make a binary corresponding to the frames he was in the bin
        tmpTimeBinary(binInds) = 1;
        
        if velFilt == 1
            tmpComboBnry = zeros(1,length(tmpTimeBinary));
            tmpComboBnry(tmpTimeBinary == 1 & velFiltBnry == 1) = 1;
            tmpTimeBinary = tmpComboBnry;
        end
        
        if sum(tmpTimeBinary) > 0 %if rat was in this spatial bin
            
            
            timeChunks = bwconncomp(tmpTimeBinary,4); %find chunks of time where rat was in this spatial bin
            %NOTE: I'm doing it this 'chunk' way vs. just multiplying # of frames by time-per-frame because it helps in finding spikes below
            %      where I search for all spikes that occurred within each time window
            
            cntr = 1;
            for c = 1:length(timeChunks.PixelIdxList)
                tmpChunkInds = [timeChunks.PixelIdxList{c}(1) timeChunks.PixelIdxList{c}(end)];
                tmpChunkTimes = coords(tmpChunkInds,1);
                binTimes{x,y}(cntr,:) =[tmpChunkTimes(1)-halfFt tmpChunkTimes(2)+halfFt]; %get start and end time for the chunk
                cntr = cntr + 1;
            end %for each chunk of time the rat was in this bin
            
            
            binTimeSum(x,y) = sum(diff(binTimes{x,y},[],2)); %amt of time rat spent in this spatial bin
            
        end %if the rat was in this spatial bin at all
        
    end %yRange
end %xRange


%% GET SPIKE COUNTS FOR EACH SPATIAL BIN
binSpkCnt = zeros(numXBins, numYBins);
for x = 1:numXBins
    for y = 1:numYBins
        if ~isempty(binTimes{x,y})
            for t = 1:size(binTimes{x,y},1)
                tmpTimeBnds = binTimes{x,y}(t,:);
                spkInds = find(spkTms>=tmpTimeBnds(1) & spkTms<=tmpTimeBnds(2));
                numSpks = length(spkInds);
                
                if velFilt == 1 %if only looking at firing rate while the rat is moving
                    if ~isempty(spkInds)
                        for si = 1:length(spkInds)
                            spkTmInd = find(coords(:,1)>=spkTms(si), 1, 'First'); % doing '>=' because of the half-frame time adjustment above
                            if velFiltBnry(spkTmInd) == 0 %if the rat wasn't moving when this spike went off...
                                numSpks = numSpks - 1; %don't count this spike
                            end
                        end
                    end
                end
                
                binSpkCnt(x,y) = binSpkCnt(x,y) + numSpks;
            end %each chunk of time
        end %if the rat was in this spatial bin
    end %y
end %x



%% FLIP binTimeSum AND binSpkCnt AROUND SO THEY'RE LINED UP WITH COORDS THE RIGHT WAY
%    This happens because coords is referenced to the bottom left corner of a plot but
%    the x & y referencing above starts in the top left (i.e., binTimeSum(1,1) = top left corner).
binTimeSum = flip(rot90(binTimeSum));
binSpkCnt = flip(rot90(binSpkCnt));




%% IF MINIMUM TIME IN SPATIAL BIN CRITERIA BEING ENFORCED...
%   THEN REMOVE SPIKES FROM ANY SPATIAL BIN THAT RAT VISITED < minBinDur ms
if durCrit == 1
    binSpkCnt(binTimeSum<(minBinDur/1000)) = 0;
end




%% CALCULATE THE RATE-MAP
rateMap = binSpkCnt ./ binTimeSum;
if plotOrNot == 1
    
    %ratemap
    subplot(1,2,2);
    plot_ratemap(rateMap, spatBinSz, 0, 1); 
    caxis([-0.01 max(rateMap(:))]); %MD: CLIM TO CAXIS
    colorbar
    xlabel('Position (cm)');
    ylabel('Position (cm)');
    title('Ratemap');
    set(gca, 'FontSize', 14)
    axis square;
    
    
    figure('Position', [411 90 1109 389])
    %spkCnts
    subplot(1,2,1);
    plot_ratemap(binSpkCnt, spatBinSz);
    colorbar
    xlabel('Position (cm)');
    ylabel('Position (cm)');
    title('Spike Counts');
    set(gca, 'FontSize', 14)
    axis square;
    
    
    %timePerBin
    subplot(1,2,2);
    plot_ratemap(binTimeSum, spatBinSz);
    colorbar
    xlabel('Position (cm)');
    ylabel('Position (cm)');
    title('Time Per Bin');
    set(gca, 'FontSize', 14)
    axis square;
    
end


%% ASSIGN OUTPUT
spkCnt = binSpkCnt;
timePerBin = binTimeSum;

end %function




