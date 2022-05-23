function group = fmr1CircTrack_x_makeStruct
curDir = pwd;
structNickName = 'nick';

% group = fmr1CircTrack_1_buildDataStruct; %calls fmr1CircTrack_0_
% group = fmr1CircTrack_2_attachPfs(group);
% group = fmr1CircTrack_3_tagLaps(group);
% group = fmr1CircTrack_4_detectSequences(group);
% group = fmr1CircTrack_5_addSleepInfoToStruct(group);
% group(2).rat.day.sleep(1) = [];
% group(2).rat.day.sleep(4) = [];

group = fmr1CircTrack_6_addReplayEventsToStruct(group);

%cd('E:\FMR1_CIRCTRACK\DATA_STRUCTS')
cd("C:\Users\nick\Projects\DATA_STRUCTS")

tmpDate = clock;
if tmpDate(2) < 10
    strDate = [num2str(tmpDate(1)) '0' num2str(tmpDate(2)) num2str(tmpDate(3))];
else
    strDate = [num2str(tmpDate(1)) num2str(tmpDate(2)) num2str(tmpDate(3))];
end %if we need to add a 0 to month

save(['dataStruct_postFxn6_' structNickName strDate], 'group')
close all
cd(curDir)

% UPDATE 7/21: Added detect sequences code
% UPDATE 8/13/21: Changed from 5 deg to 4 deg bin size
% UPDATE 9/23/21: Detect sequences/decoding with two potenital methods
% UPDATE 9/28/21: Changed min firing rate for unit to be included to 0.5 Hz
%   Added day 3 for Rat326Z
%   Fixed seq times in detect_sequence_events
% UPDATE 9/30/21: Changed min firing rate for unit bac to 1 Hz
%   Changed code to detect SWRs using Ernie method (DetectRipples_v4)
% UPDATE 10/1/21: Changed fxn 5 to only add the sleep info (spikes, coords,
%   etc). Fxn 6 now is for adding replay events via the population events
%   code and the detected SWRs from the LFP.

end %function