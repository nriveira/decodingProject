function replayEvents = nick_analyzeSync(replayEvents, saveDir)
    for g = 1:length(replayEvents)
        normDay = 1;
    
        replayEvents(g).lengthCounts = zeros(length(replayEvents(g).events),1);
        for e = 1:length(replayEvents(g).events)
            replayEvents(g).lengthCounts(e) = length(replayEvents(g).events(e).lfpSWR)/2;
        end
    end
    
    figure()
    for g = 1:length(replayEvents)
        subplot(2,1,g); axis square;
        wavelet = nan(250,replayEvents(g).max_length,length(replayEvents(g).events));
    
        % Change Normalization
        for e = 1:length(replayEvents(g).events)
            if(normDay==1)
                re = replayEvents(g).events(e).wavelet_power_dayPre;
            else
                re = replayEvents(g).events(e).wavelet_power_sleepPre;
            end
    
            sizeWP = size(re);
            wavelet(:,1:sizeWP(2),e) = re;
        end
        waveletValues = mean(wavelet, 3, 'omitnan');
       
        x = (-800:length(waveletValues(1,:))-800)./2;
        y = 1:250;
        
        imagesc(x,y,waveletValues);
        colormap('jet');
        set(gca,'YDir','normal')
    
        c = colorbar;   
        ylabel(c, 'Power (z-scored)')
    
        xline(0, 'k--')
        xlim([-400, 400])
        xlabel('Time since SWR Detection (ms)')
        ylabel('Frequency (Hz)')
        title(strcat(replayEvents(g).name, " during Sleep"))
    end
    cd(saveDir); saveas(gcf,'replayHeatmaps.png')
end %function