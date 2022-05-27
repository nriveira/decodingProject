function group = fmr1CircTrack_0_hardCodeTestData(group)
% function group = fmr1CircTrack_0_hardCodeTestData(group)
%
% Function adds the names, dates, reward locs, and theta tet to use to the
% group structure. It's just separated to this function vs it's parent because
% it takes less space.
%
% JBT
% 03/21
% Colgin Lab




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%                 KO RATS                 %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% RAT 316
group(2).rat(1).name = 'rat316';
group(2).rat(1).day(1).name = '2020-11-08'; %Day 1
group(2).rat(1).day(1).rewLocs = [90 270]; % N/S
group(2).rat(1).day(1).thetaTet = 4; %5 units, highest amp/cleanest theta
% -->>>>> Theta tetrodes chosen by manual inspection of LFPs across all tets for each day
%          - Evaluated for # of units, clarity of theta, and amplitude of theta

group(2).rat(1).day(2).name = '2020-11-09'; %Day 2
group(2).rat(1).day(2).rewLocs = [0 180]; % E/W
group(2).rat(1).day(2).thetaTet = 4; %7 units, highest amp/cleanest theta

group(2).rat(1).day(3).name = '2020-11-10'; %Day 3
group(2).rat(1).day(3).rewLocs = [135 315]; % NW/SE
group(2).rat(1).day(3).thetaTet = 4; %1 units, highest amp/cleanest theta

group(2).rat(1).day(4).name = '2020-11-11'; %Day 4
group(2).rat(1).day(4).rewLocs = [45 225]; % NE/SW
group(2).rat(1).day(4).thetaTet = 7; %1 units, highest amp/cleanest theta

group(2).rat(1).day(5).name = '2020-11-12'; %Day 5
group(2).rat(1).day(5).rewLocs = [125 305]; % NNW/SSE
group(2).rat(1).day(5).thetaTet = 7; %1 units, highest amp/cleanest theta

group(2).rat(1).day(6).name = '2020-11-16'; %Day 6
group(2).rat(1).day(6).rewLocs = [67.5 247.5]; % NNE/SSW
group(2).rat(1).day(6).thetaTet = 2; %only option

group(2).rat(1).day(7).name = '2020-11-17'; %Day 7
group(2).rat(1).day(7).rewLocs = [145 325]; % NW/SE
group(2).rat(1).day(7).thetaTet = 2; %only option

group(2).rat(1).day(8).name = '2020-11-19'; %Day 8
group(2).rat(1).day(8).rewLocs = [60 240]; % NE/SW
group(2).rat(1).day(8).thetaTet = 3; %only option

group(2).rat(1).day(9).name = '2020-11-21'; %Day 9
group(2).rat(1).day(9).rewLocs = [22.5 202.5]; % ENE/WSW
group(2).rat(1).day(9).thetaTet = 2; %only option

group(2).rat(1).day(10).name = '2020-11-23'; %Day 10
group(2).rat(1).day(10).rewLocs = [160 340]; % WNW/ESE
group(2).rat(1).day(10).thetaTet = 2; %only option

group(2).rat(1).day(11).name = '2020-11-27'; %Day 11
group(2).rat(1).day(11).rewLocs = [90 270]; %N/S
group(2).rat(1).day(11).thetaTet = 2; %2 units, cleanest (others not good looking)

group(2).rat(1).day(12).name = '2020-11-28'; %Day 12
group(2).rat(1).day(12).rewLocs = [0 180]; %E/W
group(2).rat(1).day(12).thetaTet = 2; %3 tetNums, cleanest

group(2).rat(1).day(13).name = '2020-11-30'; %Day 13
group(2).rat(1).day(13).rewLocs = [45 225]; %NE/SW
group(2).rat(1).day(13).thetaTet = 4; %5 units, cleanest, high amp

%% RAT 330
group(2).rat(2).name = 'rat330';
group(2).rat(2).day(1).name = '2021-03-04'; %Day 1
group(2).rat(2).day(1).rewLocs = 45; % NE (only 1 rew loc)
group(2).rat(2).day(1).thetaTet = 3; %8 units, tied for most units, cleanest theta

group(2).rat(2).day(2).name = '2021-03-08'; %Day 2 FOR MEG for Emma this is different
group(2).rat(2).day(2).rewLocs = [135 315]; % NW/SE
group(2).rat(2).day(2).thetaTet = 3; %selected by Emma

group(2).rat(2).day(3).name = '2021-03-15'; %Day 3 FOR MEG for Emma this is different
group(2).rat(2).day(3).rewLocs = [135 315]; % NW/SE
group(2).rat(2).day(3).thetaTet = 4; %selected by Emma




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%                  WT RATS                %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% RAT 113
% group(1).rat(1).name = 'rat113';
% group(1).rat(1).day(1).name = '2017-01-29'; %Day 1
% group(1).rat(1).day(1).rewLocs = [0 180]; % E/W
% group(1).rat(1).day(1).thetaTet = 12; %4 units (least #), but by far the cleanest theta

%% RAT 326Z

group(1).rat(1).name = 'rat326Z';
group(1).rat(1).day(1).name = '2021-06-14A'; %Day 1
group(1).rat(1).day(1).rewLocs = [0 180]; % E/W
group(1).rat(1).day(1).thetaTet = [6]; %10 units

group(1).rat(1).day(2).name = '2021-06-14B'; %Day 2
group(1).rat(1).day(2).rewLocs = [45 225]; %NE/SW AKA red and blue legs
group(1).rat(1).day(2).thetaTet = [5]; %4 units, but least amount of noise

group(1).rat(1).day(3).name = '2021-06-15'; %Day 3
group(1).rat(1).day(3).rewLocs = [135 315]; %NW/SE AKA green and yellow legs
group(1).rat(1).day(3).thetaTet = [5]; %EXPLAIN HERE









