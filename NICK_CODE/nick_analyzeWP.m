function replayEvents = nick_analyzeWP(group, replayEvents, saveDir)  
    normDay = 1;
    %% Estimating power on all replay events
    % Wavelet method
    figure(); hold on; axis square;
    colors = 'gr';
    for g = 1:length(group)
        if(normDay == 1)
            samples_wvt = zeros(length(replayEvents(g).events),length(replayEvents(g).events(1).wavelet_power_day));
            for i = 1:length(replayEvents(g).events)
                samples_wvt(i,:) = mean(replayEvents(g).events(i).wavelet_power_day,2);
            end
        else
            samples_wvt = zeros(length(replayEvents(g).events),length(replayEvents(g).events(1).wavelet_power_sleep));
            for i = 1:length(replayEvents(g).events)
                samples_wvt(i,:) = mean(replayEvents(g).events(i).wavelet_power_sleep,2);
            end
        end
        
        [ci_wvt, ~] = bootci(2000, {@mean, samples_wvt});
    
        x = 1:250;
        patch([x fliplr(x)], [ci_wvt(1,:) fliplr(ci_wvt(2,:))], colors(g),'facealpha',0.5,'LineStyle','none')
        plot(x,mean(ci_wvt),colors(g))
    end
    title('All Events (Wavelet)')
    xlim([25,250])
    xlabel('Frequency (Hz)')
    ylabel('Power (Z-score)')
    legend({'WT','','KO',''})
    cd(saveDir); saveas(gcf,'waveletAllEvents.png')
    
    %% Sorting forward vs reverse
    for g = 1:length(group)
        replayEvents(g).forward = replayEvents(g).events([replayEvents(g).events(:).slope] > 0);
        replayEvents(g).reverse = replayEvents(g).events([replayEvents(g).events(:).slope] < 0);
    end
    
    % Forward Wavelet
    figure(); hold on; axis square;
    colors = 'gr';
    for g = 1:length(group)
        if(normDay == 1)
            samples_wvt = zeros(length(replayEvents(g).forward),length(replayEvents(g).forward(1).wavelet_power_day));
            for i = 1:length(replayEvents(g).forward)
                samples_wvt(i,:) = mean(replayEvents(g).forward(i).wavelet_power_day,2);
            end
        else
            samples_wvt = zeros(length(replayEvents(g).forward),length(replayEvents(g).forward(1).wavelet_power_sleep));
            for i = 1:length(replayEvents(g).forward)
                samples_wvt(i,:) = mean(replayEvents(g).forward(i).wavelet_power_sleep,2);
            end
        end
        [ci_wvt, bootstat] = bootci(2000, {@mean, samples_wvt});
    
        x = 1:250;
        patch([x fliplr(x)], [ci_wvt(1,:) fliplr(ci_wvt(2,:))], colors(g),'facealpha',0.5,'LineStyle','none')
        plot(1:250,mean(ci_wvt),colors(g))
    end
    
    title('Forward (Wavelet)')
    xlim([25,250])
    xlabel('Frequency (Hz)')
    ylabel('Power (Z-score)')
    legend({'WT','','KO',''})
    (saveDir); saveas(gcf,'waveletForwardEvents.png')
    
    % Reverse wavelet
    figure(); hold on; axis square;
    colors = 'gr';
    for g = 1:length(group)
        if(normDay == 1)
            samples_wvt = zeros(length(replayEvents(g).reverse),length(replayEvents(g).reverse(1).wavelet_power_day));
            
            for i = 1:length(replayEvents(g).reverse)
                samples_wvt(i,:) = mean(replayEvents(g).reverse(i).wavelet_power_day,2);
            end
        else
            samples_wvt = zeros(length(replayEvents(g).reverse),length(replayEvents(g).reverse(1).wavelet_power_sleep));
            
            for i = 1:length(replayEvents(g).reverse)
                samples_wvt(i,:) = mean(replayEvents(g).reverse(i).wavelet_power_sleep,2);
            end 
        end
        
        [ci_wvt, ~] = bootci(2000, {@mean, samples_wvt});
    
        x = 1:250;
        patch([x fliplr(x)], [ci_wvt(1,:) fliplr(ci_wvt(2,:))], colors(g),'facealpha',0.5,'LineStyle','none')
        plot(x,mean(ci_wvt),colors(g))
    end
    
    title('Reverse (Wavelet)')
    xlim([25,250])
    xlabel('Frequency (Hz)')
    ylabel('Power (Z-score)')
    legend({'WT','','KO',''})
    cd(saveDir); saveas(gcf,'waveletReverseEvents.png')
end