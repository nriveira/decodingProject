function runEvents = nick_analyzeRun(runEvents, saveDir)
    for g = 1:length(runEvents)
        normDay = 1;
    
        runEvents(g).lengthCounts = zeros(length(runEvents(g).events),1);
        for e = 1:length(runEvents(g).events)
            runEvents(g).lengthCounts(e) = length(runEvents(g).events(e).cutTheta);
        end
        runEvents(g).max_length = max(runEvents(g).lengthCounts);
    end
    
    figure()
    for g = 1:length(runEvents)
        subplot(2,1,g); axis square;
        wavelet = nan(250,runEvents(g).max_length, length(runEvents(g).events));
    
        % Change Normalization
        for e = 1:length(runEvents(g).events)
            if(normDay==1)
                re = runEvents(g).events(e).wavelet_power_day;
            else
                re = runEvents(g).events(e).wavelet_power_begin;
            end
    
            sizeWP = size(re);
            wavelet(:,1:sizeWP(2),e) = re;
        end
        waveletValues = mean(wavelet, 3, 'omitnan');
        % X is degrees in Theta Cycle?
        x = (1:size(waveletValues,2));
        y = 1:250;
        
        imagesc(x,y,waveletValues);
        colormap('jet');
        set(gca,'YDir','normal')  
        xlabel('Theta (Degrees)')
        ylabel('Frequency (Hz)')
        title(strcat(runEvents(g).name,  " during Run"))

        c = colorbar;
        ylabel(c, 'Power (z-scored)')
    end
    cd(saveDir); saveas(gcf,'runHeatmaps.png')
end %function