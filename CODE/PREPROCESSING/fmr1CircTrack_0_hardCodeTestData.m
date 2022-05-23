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

%% RAT 330
group(2).rat(1).name = 'rat330';
group(2).rat(1).day(1).name = '2021-03-04'; %Day 1
group(2).rat(1).day(1).rewLocs = 45; % NE (only 1 rew loc)
group(2).rat(1).day(1).thetaTet = 3; %8 units, tied for most units, cleanest theta

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                         %
%                  WT RATS                %
%                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% RAT 326Z
group(1).rat(1).name = 'rat326Z';
group(1).rat(1).day(1).name = '2021-06-14B'; %Day 2
group(1).rat(1).day(1).rewLocs = [45 225]; %NE/SW
group(1).rat(1).day(1).thetaTet = [5]; %4 units, but least amount of noise









