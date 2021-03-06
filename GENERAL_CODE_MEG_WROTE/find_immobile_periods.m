function [immTimes, timeImm] = find_immobile_periods(coords)
% function [immTimes, timeImm] = find_immobile_periods(coords)
% 
%   PURPOSE:
%       Identify the periods of time where the animal can be considered
%       immobile. This is when the animal is moving > 5 cm/s for a minimum
%       of 5 seconds. 
% 
%   INPUT:
%       coords = struct from read_in_coords
% 
%   OPTIONS:
%       Internally, you can change the max speed and minimum duration for
%       the animal to be considered immobile. The runspeed is also
%       smoothed, and you can change the standard deviation for smoothing
%       (default is 250/6 ms, with 250 ms as the window).
% 
%   OUTPUT:
%       immTimes = struct with times when the animal is immobile
%           (:,1) = start of immobile periods
%           (:,2) = end of immobile periods
%       timeImm = length of time (in s) animal was immobile
%
% MM Donahue
% 04/2020
% Colgin Lab

%% SET PARAMETERS/INITIALIZE

maxSpd = 5; %cm/s
minDur = 1; %seconds of immobility to be included

sampRate = 29.97; %frames/sec for camera
minSamp = minDur * sampRate;

gWinDur = 250; %ms, for smoothing run speed
gWinStd = 250/6; %ms, for smoothing run speed

immTimes = []; %initialize

%% GET TIMES

timeImm = 0;

instRunSpd = get_runspeed(coords);
smRunSpd = smooth_runspeed(instRunSpd, gWinDur, gWinStd);

slowInds = smRunSpd(:,2) < maxSpd;

chunks = bwconncomp(slowInds,4); %find where multiple consecutive frames of slow speed are

for ch = 1:chunks.NumObjects %for all potenital immobile periods
    curChunk = chunks.PixelIdxList{ch};
    if length(curChunk) > minSamp %see if it is at least minDur seconds
       
       startInd = curChunk(1,1);
       endInd = curChunk(end,1);
       
       startTm = coords(startInd,1); %convert to seconds
       endTm = coords(endInd,1);
       
       immTimes = [immTimes; startTm endTm]; %save it for output
       timeImm = timeImm + (endTm - startTm);

    end %if long enough
end %all detected conn comps



end %function