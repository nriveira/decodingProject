function fmr1CircTrack_x_plot_spikeTimeAutocorrelations(group)
% function fmr1CircTrack_x_plot_spikeTimeAutocorrelations(group)
% 
% PURPOSE:
%   Plot spike time autocorrelations for cells from WT and KO rats.
% 
% INPUTS:
%   group struct
% 
% OUTPUTS:
%   Figures.
% 
% OPTIONS:
%   See code for options to plot by individual rats and days. Options to
%   save figs.
% 
% MMD
% 07/2021
% Colgin Lab


%% OPTIONS

plotByDay = 0; %1 to plot mean for every day, 0 to not
plotByRat = 0; %1 to plot mean for every rat, 0 to not

saveOrNot = 0; %1 to save, 0 to not

saveOrNot = 1;
saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\spikeTimeAutocorrelations';
curDir = pwd;

cols = {'Blue', 'Red'};

%% INITIALIZE

totWin = 2000; %ms
binSz = 10/1000; %10 ms

autoCorr = cell(2,1);

%% GET DATA

for g = 1:2
    fprintf('%s\n', group(g).name)
    for r = 1:length(group(g).rat)
        fprintf('\t%s\n', group(g).rat(r).name)
        autoCorrByRat = [];
        
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day))
            
            autoCorrByDay = [];
            
            for u = 1:length(group(g).rat(r).day(d).xAllBeginUnitInfo)
                fprintf('\t\t\tUnit %d/%d\n', u, length(group(g).rat(r).day(d).xAllBeginUnitInfo))
                autoCorrForUnit = [];
                
                for b = 1:4
                    
                    spkTms = group(g).rat(r).day(d).begin(b).unit(u).spkTms;
                    
                    timeStart = spkTms(1);
                    timeEnd = spkTms(end);
                    totTime = timeEnd - timeStart; %time in seconds
                    if totTime ~= 0
                        
                        numBins = round(totTime/binSz);
                        
                        tmpBinSpks = histcounts(spkTms, timeStart:binSz:timeEnd);
                        keyboard
                        try
                            tmpAutoCorr = xcorr(tmpBinSpks, 100, 'coeff');
                        catch
                            keyboard
                        end
                        autoCorrForUnit = [autoCorrForUnit; tmpAutoCorr];
                    end
                end %begin
                
                autoCorrByRat = [autoCorrByRat; mean(autoCorrForUnit,1)];
                autoCorrByDay = [autoCorrByDay; mean(autoCorrForUnit,1)];
                autoCorr{g} = [autoCorr{g}; mean(autoCorrForUnit,1)];
                
            end %unit
            
            %Make within day figure
            if plotByDay == 1
                figName = [group(g).name '_' group(g).rat(r).name '_' group(g).rat(r).day(d).name];
                figure('Name', figName)
                
                %             meanSpksPerBin = mean(autoCorrByRat,1);
                %             convFiringRate = meanSpksPerBin * (1000/binSz);
                %             plot(convFiringRate)
                plot(mean(autoCorrByDay,1), 'Color', rgb(cols{g}), 'LineWidth', 1.5)
                
                xlim([0 200])
                xticks([0:50:200])
                xticklabels({'-1000' '-500' '0' '500' '1000'})
                xlabel('Time (ms)')
                
                ttl = {[group(g).name ' - ' group(g).rat(r).name ' - ' group(g).rat(r).day(d).name]; ['n = ' num2str(size(autoCorrByDay,1)) ' cells']};
                title(ttl)
                
                ylim([0 .15])
                ylabel('Correlation coefficient')
                yRange = get(gca, 'YLim');
                han = line([100 100], yRange);
                set(han, 'LineStyle', '--', 'Color', [0 0 0]);
                
                if saveOrNot == 1
                    cd(saveDir)
                    saveas(gcf, figName, 'epsc')
                    saveas(gcf, figName, 'png')
                    saveas(gcf, figName, 'fig')
                    cd(curDir)
                end %save option
            end %if plot by day
            
        end %day
        
        % Make all days for this rat figure
        if plotByRat == 1
            figName = [group(g).name '_' group(g).rat(r).name '_allDays'];
            figure('Name', figName)
            
            plot(mean(autoCorrByRat,1), 'Color', rgb(cols{g}), 'LineWidth', 1.5)
            
            
            xlim([0 200])
            xticks([0:50:200])
            xticklabels({'-1000' '-500' '0' '500' '1000'})
            xlabel('Time (ms)')
            
            ttl = {[group(g).name ' - ' group(g).rat(r).name]; ['n = ' num2str(size(autoCorrByRat,1)) ' cells']};
            title(ttl)
            
            ylim([0 0.15])
            ylabel('Correlation coefficient')
            yRange = get(gca, 'YLim');
            han = line([100 100], yRange);
            set(han, 'LineStyle', '--', 'Color', [0 0 0]);
            
            if saveOrNot == 1
                cd(saveDir)
                saveas(gcf, figName, 'epsc')
                saveas(gcf, figName, 'png')
                saveas(gcf, figName, 'fig')
                cd(curDir)
            end %save option
        end %plot by rat
        
    end %rat
end %group

%% FIG 1(ish) - 


figtitle = 'SpikeTimeAutoCorrelations';
figure('Name', figtitle, 'Position', [576 528 861 420])


for g = 1:2
    subplot(1,2,g)
    
    plot(mean(autoCorr{g},1), 'Color', rgb(cols{g}), 'LineWidth', 1.5)
    
    
    xlim([0 200])
    xticks([0:50:200])
    xticklabels({'-1000' '-500' '0' '500' '1000'})
    xlabel('Time (ms)')
    
    ttl = {[group(g).name]; ['n = ' num2str(size(autoCorr{g},1)) ' cells']};
    title(ttl)
    
    ylim([0 0.1])
    ylabel('Correlation coefficient')
    yRange = get(gca, 'YLim');
    han = line([100 100], yRange);
    set(han, 'LineStyle', '--', 'Color', [0 0 0]);
    
    
end %group

if saveOrNot == 1
    cd(saveDir)
    saveas(gcf, figtitle, 'epsc')
    saveas(gcf, figtitle, 'png')
    saveas(gcf, figtitle, 'fig')
    cd(curDir)
end %save option


keyboard
end %function