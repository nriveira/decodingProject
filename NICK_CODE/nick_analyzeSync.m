load("C:\Users\nick\Projects\DATA_STRUCTS\replayEvents_postFxn6_nick20220413.mat")

for g = 1:length(replayEvents)
    replayEvents(g).lengthCounts = zeros(length(replayEvents(g).events),1);
    for e = 1:length(replayEvents(g).events)
        replayEvents(g).lengthCounts(e) = length(replayEvents(g).events(e).lfpSWR)/2;
    end
end

figure()
for g = 1:length(replayEvents)
    subplot(2,1,g); axis square;
    wavelet = nan(226,replayEvents(g).max_length+1600,length(replayEvents(g).events));
    for e = 1:length(replayEvents(g).events)
        sizeWP = size(replayEvents(g).events(e).wavelet_powerPre);
        wavelet(:,1:sizeWP(2),e) = replayEvents(g).events(e).wavelet_powerPre;
    end
    waveletValues = mean(wavelet, 3, 'omitnan');
   
    x = (-800:length(waveletValues(1,:))-800)./2;
    y = 25:250;
    
    imagesc(x,y,waveletValues);
    colormap('jet');
    set(gca,'YDir','normal')

    c = colorbar;   
    caxis([20,35])
    ylabel(c, 'Power (dB)')

    xline(0, 'k--')
    xlim([-400, 400])
    xlabel('Time since SWR Detection (ms)')
    ylabel('Frequency (Hz)')
    title(replayEvents(g).name)
end
saveas(gcf,'figures/replayHeatmaps.png')

figure()
x = [replayEvents(1).lengthCounts; replayEvents(2).lengthCounts];
g = [repmat({'WT'},202,1); repmat({'KO'},115,1)];
boxplot(x, g)
title('Event Durations')
ylabel('Time (ms)')
saveas(gcf,'figures/replayEventCounts.png')