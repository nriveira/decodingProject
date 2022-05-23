function fmr1CircTrack_x_plotPropFiringUnits(group)
% function fmr1CircTrack_x_plotPropFiringUnits(group)
%
% PURPOSE:
%  Function makes a line-graph (line for each group) showing the proportion of cells firing > 0.5 Hz
%  by location, similar to Zheng et al., BioRxiv, Figure 1b
%
% INPUT:
%  group = project uber data structure
%
% OUTPUT:
%  The figure described above
%
% JB Trimper
% 1/6/20
% Colgin Lab


binSize = 4; %degrees
binCtrs = linspace(binSize/2, 360-binSize/2, 360/binSize);
newRewLoc = 90;
[~,newRewInd] = min(abs(circ_dist(deg2rad(binCtrs), deg2rad(newRewLoc))-0));

Ns = [0 0];

for g = 1:2
    fprintf('Group %d\n', g);
    
    unitsAboveThreshBnry = [];
    
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);
        for d = 1:length(group(g).rat(r).day)
            
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day));
            if length(group(g).rat(r).day(d).xAllBeginUnitInfo) >=20
            % Get and smooth ratemaps
            rateMaps = zeros(length(group(g).rat(r).day(d).xAllBeginUnitInfo), length(group(g).rat(r).day(d).xAllBeginUnitInfo(1).smRateMap));
            binCtrs = group(g).rat(r).day(d).binCtrs; %will potentially repeat, but it doesn't matter
            for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                smRm = group(g).rat(r).day(d).xAllBeginUnitInfo(u).smRateMap;
                rateMaps(u,:) = smRm;
            end
            
            
            % Shift the ratemaps so reward locations are always on top of one another when plotting
            rewLoc = group(g).rat(r).day(d).rewLocs(1);
            [~,rewInd] = min(abs(circ_dist(deg2rad(binCtrs), deg2rad(rewLoc))));
            shiftVal = newRewInd - rewInd;
            rateMaps = circshift(rateMaps, [0 shiftVal]);
            
            % Find the number of units above threshold and concatenate that into the bigger matrix
            tmpUnitsAboveThreshBnry = zeros(size(rateMaps));
            tmpUnitsAboveThreshBnry(rateMaps>=0.5) = 1;
            if r == 1 && d == 1
                unitsAboveThreshBnry = tmpUnitsAboveThreshBnry; %if its the first rat/day for the group, declare this matrix
            else
                unitsAboveThreshBnry = [tmpUnitsAboveThreshBnry; unitsAboveThreshBnry]; %#ok -- if not, concatenate onto it
            end
            end
        end %day
    end %rat
    
    % Calculate proportions of cells firing (if there were cells)
    if ~isempty(unitsAboveThreshBnry)
        propFiring(g,:) = sum(unitsAboveThreshBnry,1) / size(unitsAboveThreshBnry,1); %#ok
        Ns(g) = size(unitsAboveThreshBnry,1);
    else
        propFiring(g,:) = zeros(1,360/binSize); %#ok
    end
    
end %group


%% PLOTTING

figure('Position', [631   244   617   531]);
hold on;

% Format, label
axis square;
ylabel('Unit Proportion');
xlabel('Aligned Angular Position (degrees)');
title('Proportion of Total Units >=0.5 Hz by position');
fix_font;
xlim([0 360]);

%Mark reward locations with gray lines
yMax = 0.75;
ylim([0 yMax]); % Markers need to go from bottom to top of the panel so you need yMax established
ln = line([newRewLoc newRewLoc], [0 yMax]);
set(ln, 'Color', [.5 .5 .5], 'LineWidth', 5);
ln = line([newRewLoc+180 newRewLoc+180], [0 yMax]);
set(ln, 'Color', [.5 .5 .5], 'LineWidth', 5);

% Plot the actual proportions by location
for g = 1:2
    ln(g+1) = plot(binCtrs, propFiring(g,:), 'LineWidth', 2);
end

legend(ln, {'Rew Locs', ['WT (n = ' num2str(Ns(1)) ')'], ['KO (n = ' num2str(Ns(2)) ')' ]}, 'Location', 'SouthEast');


end %fnctn