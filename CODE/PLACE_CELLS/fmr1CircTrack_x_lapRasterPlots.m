function fmr1CircTrack_x_lapRasterPlots(group)
% function fmr1CircTrack_x_lapRasterPlots(group)
% 
% PURPOSE:
%   Plot example spike raster plots for successive laps across all four
%   begins. Makes plots for all units. Saves in folders for each day for
%   each rat.
% 
% INPUT:
%   group struct
% 
% OUTPUT:
%   Figures that can optionally be saved (see function). 
%   Function closes figures internally, so if you want to see them comment 
%       out line 102 or add a keyboard. Makes an individual plot for each
%       unit, so be aware it could be a lot of figures.
% 
% MMD
% 7/2021
% Colgin Lab

%% OPTIONS

saveOrNot = 1;
saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\PLACE_CELLS\lapRasterPlots';

cols = {'Blue', 'Red'};

%% INITIALIZE

binSize = 5; %degrees
binCtrs = linspace(binSize/2, 360-binSize/2, 360/binSize);

cd(saveDir)

%% GET DATA/MAKE PLOTS

for g = 1:2
    cd(group(g).name)
    for r = 1:length(group(g).rat)
        if ~isfolder(group(g).rat(r).name)
            mkdir(group(g).rat(r).name)
        end %if is folder
        cd(group(g).rat(r).name)
        for d = 1:length(group(g).rat(r).day)
            
            if ~isfolder(group(g).rat(r).day(d).name)
                mkdir(group(g).rat(r).day(d).name)
            end %if is folder
            
            cd(group(g).rat(r).day(d).name)
            for u = 1:length(group(g).rat(r).day(d).xBeginUnitInfo)
                tmpID = group(g).rat(r).day(d).xBeginUnitInfo(u).ID;
                figtitle = ['TT' num2str(tmpID(1)) '_' num2str(tmpID(2))];
                
                figure('Name', figtitle)
                
                yVal = 1;
                lpCntr = 0;
                
                for b = 1:4
                    
                    spkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                    radPos = group(g).rat(r).day(d).begin(b).radPos;
                    
                    for lp = 1:size(group(g).rat(r).day(d).begin(b).lapTms,1)
                        lapStart = group(g).rat(r).day(d).begin(b).lapTms(lp,1);
                        lapEnd = group(g).rat(r).day(d).begin(b).lapTms(lp,2);
                        lapSpks = spkTms(spkTms >= lapStart & spkTms <= lapEnd);
                        
                        lapSpkPos = [];
                        
                        for st = 1:length(lapSpks)
                            spkTm = lapSpks(st);
                            
                            posInd = find(radPos(:,1) <= spkTm, 1, 'Last');
                            
                            lapSpkPos = [lapSpkPos radPos(posInd,2)];
                            line([radPos(posInd,2) radPos(posInd,2)], [yVal-.3 yVal+.3], 'Color', rgb(cols{g}))
                            
                        end %spktms
                        
                        yVal = yVal + 1;
                    end %lap
                    
                    ln = line([binCtrs(1) binCtrs(end)], [yVal-.5 yVal-.5]);
                    set(ln, 'Color', [.4 .4 .4], 'LineStyle', '--');
                    
                    text(binCtrs(end)+5, (size(group(g).rat(r).day(d).begin(b).lapTms,1)/2)+lpCntr, ['Begin ' num2str(b)]);
                    lpCntr = lpCntr + size(group(g).rat(r).day(d).begin(b).lapTms,1);
                end %begin
                
                xlim([0 360])
                xlabel('Position (deg)')
                
                ylim([0.45 yVal+0.45])
                ylabel('Lap Number')
                
                title(['Unit: ' 'TT' num2str(tmpID(1)) '\_' num2str(tmpID(2))])
                
                if saveOrNot == 1
                    saveas(gcf, figtitle, 'epsc')
                    saveas(gcf, figtitle, 'png')
                    saveas(gcf, figtitle, 'fig')
                end %save option
                
            end %unit
            close all
            cd ../
            
        end %day
        cd ../
    end %rat
    cd ../
end %group


end %function