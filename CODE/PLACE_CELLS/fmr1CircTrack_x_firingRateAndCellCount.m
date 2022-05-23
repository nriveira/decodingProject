function fmr1CircTrack_x_firingRateAndCellCount(group)
% function fmr1CircTrack_x_firingRateAndCellCount(group)
%
% PURPOSE:
%  To plot firing rate histograms for each group/day/rat with unit counts
%  in the subplot titles. 
%
% INPUT:
%  group = project uber data struct
%
% OUTPUT:
%   Plot described under PURPOSE, auto-saved in 'saveDir' at top of function if desired
%
% JB Trimper
% 01/2021
% Colgin Lab

savePlots = 0;
saveDir = 'C:\Users\John\Desktop\LAB_STUFF\PROJECTS\FMR1_CircTrack\RESULTS\firingRates';

for g = 1:2
    fprintf('Group %d\n', g);
    
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);
        
        figure('Position', [205.8 183.4 859.2 594.4], 'Name', [group(g).name ': Rat ' num2str(r)]);
        
        numU = zeros(1,length(group(g).rat(r).day));
        
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day));
            
            numU(d) = length(group(g).rat(r).day(d).xAllBeginUnitInfo);
            spkCnts = zeros(1, length(group(g).rat(r).day(d).xAllBeginUnitInfo));
            
            sesDur = 0;
            for b = 1:4
                fprintf('\t\t\tBegin %d\n', b)
                
                startBegin = group(g).rat(r).day(d).begin(b).coords(1,1);
                endBegin = group(g).rat(r).day(d).begin(b).coords(end,1);
                sesDur = sesDur + endBegin - startBegin;
                
                for u = 1:length(group(g).rat(r).day(d).begin(b).unit)
                    spkCnts(u) = spkCnts(u) + length(group(g).rat(r).day(d).begin(b).unit(u).spkTms);
                end
                
            end %begin
            
            firingRates = spkCnts / sesDur;
            
            subplot(4,5,d);
            histogram(firingRates, 'BinWidth', 0.25);
            title(['Day' num2str(d) ' (n = ' num2str(numU(d)) ')']);
            if mod(d,5) == 1
                ylabel('Unit Count');
            end
            if d>length(group(g).rat(r).day)-5
                xlabel('Firing Rate (Hz)');
            end
            yBnds = get(gca, 'Ylim');
            yMax = yBnds(2)+.2*yBnds(2);
            ylim([0 yMax]);
            xBnds = get(gca, 'XLim');
            if xBnds(2) <= 2
                xBnds(2) = 2;
            end
            xlim([-0.05 xBnds(2)]);
            
        end %day
        
        if savePlots == 1
            curDir = pwd;
            cd(saveDir);
            print([group(g).name '_' group(g).rat(r).name], '-dpng');
            cd(curDir);
        end
        
    end %rat
end %group


end %fnctn






