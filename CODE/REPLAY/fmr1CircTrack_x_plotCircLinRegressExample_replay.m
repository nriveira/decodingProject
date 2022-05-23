function fmr1CircTrack_x_plotCircLinRegressExample_replay(group)
% function fmr1CircTrack_x_plotCircLinRegressExample_replay(group)

%% OPTIONS
figsPerDay = 20;

saveOrNot = 1; %to save figs

r2Thresh = 0.8; %set it a little higher, because this is just for plotting examples

%% INITIALIZE

curDir = pwd; %for returning after saving figs
saveDir = 'E:\FMR1_CIRCTRACK\RESULTS\REPLAY\circLinReg_examplePlots';
cd(saveDir)

degBinCtrs = 2:4:360;
radBinCtrs = deg2rad(degBinCtrs);

bound = 20; %for replay event

bayesStep = 0.01;

%% GET DATA

for g = 1:2
    for r = 1:length(group(g).rat)
        for d = 1:length(group(g).rat(r).day)
            figTtlBase = [group(g).name '_' group(g).rat(r).name '_' group(g).rat(r).day(d).name ];
            figCntr = 0; %reset for this day
            
            for s = 2:length(group(g).rat(r).day(d).sleep) %after experience, so just sleep 2 on
                
                if isempty(group(g).rat(r).day(d).sleep(s).unit) %if there is anything in the sleep
                    continue
                end %anything in sleep
                popEvents = group(g).rat(r).day(d).sleep(s).popEv; %shorten variable name
                
                for i = 1:length(popEvents)
                    
                    if popEvents(i).r2 < r2Thresh || isnan(popEvents(i).r2)
                        continue
                    end %doesn't meet threshold
                    figCntr = figCntr + 1;
                    if figCntr > figsPerDay
                        break
                    end %made enough figs for today
                    
                    figtitle = [figTtlBase '_' num2str(figCntr)];
                    figure('Name', figtitle, 'Position', [343 387 1159 591])
                    
                    pxn = popEvents(i).pxn;
                    
                    subplot(1,2,1)
                    imagesc(0:0.01:0.01*size(pxn,2), 1:4:360, pxn)
                    axis xy
                    axis square
                    colormap(hot)
                    ylabel('Position (?)')
                    xlabel('Time (s)')
                    
                    subplot(1,2,2)
                    lh = zeros(1,2);
                    [maxVals, maxInds ] = max(pxn);
                    radVals = radBinCtrs(maxInds);
                    
                    binsToUse = find(~isnan(maxVals));
                    timeAx = bayesStep/2:bayesStep:bayesStep*size(pxn,2);
                     lh(1) = scatter(timeAx(binsToUse), degBinCtrs(maxInds(binsToUse)), 'MarkerEdgeColor', rgb('Black') , 'MarkerFaceColor', rgb('Black'));
                    
                    para = circ_lin_regress(timeAx(binsToUse), radVals(binsToUse), bound);
                    
                    calphase =rad2deg(2*pi*para(1,1)*timeAx + para(1,2));
                    slope = rad2deg(para(1,1)*2*pi);
                    
                    hold on;
                    lh(2) = plot(timeAx, calphase, 'b--');
                     plot(timeAx, calphase-360, 'b--')
                    plot(timeAx, calphase+360, 'b--')
                    plot(timeAx, calphase+720, 'b--')
                    
                    ylim([0 360])
                    ylabel('Position (?)')
                    axis square
                    
                    xlabel('Time (s)')
                    xlim([0 bayesStep*size(pxn,2)])
                    
                    leg = legend(lh, {'Max probable position', 'Regression line'});
                    leg.Position = [0.82 0.855 0.1458 0.0618];
                    title({['r^2 = ' num2str(round(popEvents(i).r2,2))], ['slope = ' num2str(round(slope,2)) ' ?/s']})
                    
                    if saveOrNot == 1
                        saveas(gcf, figtitle, 'epsc');
                        saveas(gcf, figtitle, 'fig');
                        saveas(gcf, figtitle, 'png');
                    end %save or not
                end %i - pop events
                
            end %sleep
            close all
        end %day
    end %rat
end %group

cd(curDir)

end %function