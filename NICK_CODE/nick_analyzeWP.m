load("C:\Users\nick\Projects\DATA_STRUCTS\dataStruct_postFxn6_nick20220411.mat")

for g = 1:length(group)
    max_length = 0;
    min_length = 1000;

    fprintf('Group %d\n', g);
    replays(g).name = group(g).name;
    replays(g).events = [];
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day));
            for s = 1:length(group(g).rat(r).day(d).sleep)

                replayInds = find([group(g).rat(r).day(d).sleep(s).popEv(:).r2] >= 0.5);

                for ri = 1:length(replayInds)
                    replays(g).events = [replays(g).events; group(g).rat(r).day(d).sleep(s).popEv(replayInds(ri))];
                    len = length(group(g).rat(r).day(d).sleep(s).popEv(replayInds(ri)).lfps);

                    if(len<min_length)
                        min_length = len;
                    elseif(len > max_length)
                        max_length = len;
                    end

                end
            end
        end
    end

    replays(g).max_length = max_length;
    replays(g).min_length = min_length;
end

%% Estimating power on all replay events
% Wavelet method
figure(); hold on; axis square;
colors = 'gr';
for g = 1:length(group)
    samples_wvt = zeros(length(replays(g).events),length(replays(g).events(1).wavelet_power));
    
    for i = 1:length(replays(g).events)
        samples_wvt(i,:) = mean(replays(g).events(i).wavelet_power,2);
    end
    
    [ci_wvt, ~] = bootci(2000, {@mean, samples_wvt});

    x = 25:250;
    patch([x fliplr(x)], [ci_wvt(1,:) fliplr(ci_wvt(2,:))], colors(g),'facealpha',0.5,'LineStyle','none')
    plot(x,mean(ci_wvt),colors(g))
end
title('All Events (Wavelet)')
xlim([25,250])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
legend({'WT','','KO',''})
saveas(gcf,'figures/waveletAllEvents.png')

% Welchs method
figure(); hold on; axis square;
colors = 'gr';
for g = 1:length(group)
    samples_wch = zeros(length(replays(g).events),length(replays(g).events(1).welch));
    
    for i = 1:length(replays(g).events)
        samples_wch(i,:) = 20*log10(mean(replays(g).events(i).welch,2));
    end
    
    [ci_wch, ~] = bootci(2000, {@mean, samples_wch});

    x = replays(g).events(1).freq';

    patch([x fliplr(x)], [ci_wch(1,:) fliplr(ci_wch(2,:))], colors(g),'facealpha',0.5,'LineStyle','none')
    plot(x,mean(ci_wch),colors(g))
end
title('All Events (Welchs Method)')
xlim([25,250])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
legend({'WT','','KO',''})
saveas(gcf,'figures/welchAllEvents.png')

%% Sorting forward vs reverse
for g = 1:length(group)
    replays(g).forward = replays(g).events([replays(g).events(:).slope] > 0);
    replays(g).reverse = replays(g).events([replays(g).events(:).slope] < 0);
end

% Forward Wavelet
figure(); hold on; axis square;
colors = 'gr';
for g = 1:length(group)
    samples_wvt = zeros(length(replays(g).forward),length(replays(g).forward(1).wavelet_power));
    
    for i = 1:length(replays(g).forward)
        samples_wvt(i,:) = 20*log10(mean(replays(g).forward(i).wavelet_power,2));
    end
    
    [ci_wvt, bootstat] = bootci(2000, {@mean, samples_wvt});

    x = 25:250;
    patch([x fliplr(x)], [ci_wvt(1,:) fliplr(ci_wvt(2,:))], colors(g),'facealpha',0.5,'LineStyle','none')
    plot(25:250,mean(ci_wvt),colors(g))
end

title('Forward (Wavelet)')
xlim([25,250])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
legend({'WT','','KO',''})
saveas(gcf,'figures/waveletForwardEvents.png')

% Forward welch
colors = 'gr';
figure(); hold on; axis square;
for g = 1:length(group)
    samples_wch = zeros(length(replays(g).forward),length(replays(g).forward(1).welch));
    
    for i = 1:length(replays(g).forward)
        samples_wch(i,:) = 20*log10(mean(replays(g).forward(i).welch,2));
    end
    
    [ci_wch, bootstat] = bootci(2000, {@mean, samples_wch});

    x = replays(g).events(1).freq';
    patch([x fliplr(x)], [ci_wch(1,:) fliplr(ci_wch(2,:))], colors(g),'facealpha',0.5,'LineStyle','none')
    plot(x,mean(ci_wch),colors(g))
end

title('Forward (Welchs Method)')
xlim([25,250])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
legend({'WT','','KO',''})
saveas(gcf,'figures/welchForwardEvents.png')

% Reverse wavelet
figure(); hold on; axis square;
colors = 'gr';
for g = 1:length(group)
    samples_wvt = zeros(length(replays(g).reverse),length(replays(g).reverse(1).wavelet_power));
    
    for i = 1:length(replays(g).reverse)
        samples_wvt(i,:) = 20*log10(mean(replays(g).reverse(i).wavelet_power,2));
    end
    
    [ci_wvt, bootstat] = bootci(2000, {@mean, samples_wvt});

    x = 25:250;
    patch([x fliplr(x)], [ci_wvt(1,:) fliplr(ci_wvt(2,:))], colors(g),'facealpha',0.5,'LineStyle','none')
    plot(25:250,mean(ci_wvt),colors(g))
end

title('Reverse (Wavelet)')
xlim([25,250])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
legend({'WT','','KO',''})
saveas(gcf,'figures/waveletReverseEvents.png')

% Reverse welch
figure(); hold on; axis square;
colors = 'gr';
for g = 1:length(group)
    samples_wch = zeros(length(replays(g).reverse),length(replays(g).reverse(1).welch));
    
    for i = 1:length(replays(g).reverse)
        samples_wch(i,:) = 20*log10(mean(replays(g).reverse(i).welch,2));
    end
    
    [ci_wch, bootstat] = bootci(2000, {@mean, samples_wch});

    x = replays(g).events(1).freq';
    patch([x fliplr(x)], [ci_wch(1,:) fliplr(ci_wch(2,:))], colors(g),'facealpha',0.5,'LineStyle','none')
    plot(x,mean(ci_wch),colors(g))
end

title('Reverse (Welchs Method)')
xlim([25,250])
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
legend({'WT','','KO',''})
saveas(gcf,'figures/welchReverseEvents.png')