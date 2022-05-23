function fmr1CircTrack_x_makeReplayEventPlots(allData, dataBySlope)


cols = {'Blue', 'Red'};
group(1).name = 'WT';
group(2).name = 'KO';

subplot(2,2,1)
jitter = 0.2;
h = dotplot(1:2, allData, jitter, [rgb(cols{1}); rgb(cols{2})], [rgb('Black'); rgb('Black')]);

meanData = cellfun(@mean, allData);
labs = cell(1,2);
for g = 1:2
    line([g-0.25 g+0.25], [meanData(g) meanData(g)], 'Color', 'Black', 'LineWidth', 3) %mean lines should be longer than the jitter so you can see them
    
    labs{g} = [group(g).name ' n = ' num2str(length(allData{g})) ' events'];
end

xticklabels(labs)

subplot(2,2,3)
bgraph = bar(meanData, 'FaceColor', 'Flat');

for g = 1:2 %both types
    bgraph.CData(g,:) = rgb(cols{g});
end %types

semData = cellfun(@semfunct, allData);
hold on;
er = errorbar(1:2, meanData, semData, semData);
er.Color = [0 0 0];
er.LineStyle = 'none';

xticklabels({group(1).name, group(2).name})

subplot(2,2,2)
xMap = [1.25 3.25; 1.75 3.75];

meanDatabySlope = cellfun(@mean, dataBySlope);
for sInd = 1:2
    jitter = 0.1;
    
    h = dotplot(xMap(:,sInd), dataBySlope(:,sInd), jitter, [rgb(cols{1}); rgb(cols{2})], [rgb('Black'); rgb('Black')]);
    
    for g = 1:2
        line([xMap(g,sInd)-0.2 xMap(g,sInd)+0.2], [meanDatabySlope(g,sInd) meanDatabySlope(g,sInd)], 'Color', 'Black', 'LineWidth', 3);
    end %group
    
end %slope directinos

xlim([0.75 4.25])
xticks([mean(xMap(:,1)) mean(xMap(:,2))])
xticklabels({'Forward', 'Reverse'})
xlabel('Replay event slope direction')


legend({group(1).name, group(2).name}, 'Location', 'northeastoutside')

subplot(2,2,4)

meanDatabySlope = meanDatabySlope';
semDatabySlope = cellfun(@semfunct, dataBySlope)';

bgraph = bar(meanDatabySlope, 'FaceColor', 'Flat');
hold on;

ngroups = 2;
nbars = 2;
groupwidth = min(0.8, nbars/(nbars + 1.5));
for g = 1:2
    bgraph(g).CData = rgb(cols{g});
    
    x = (1:ngroups) - groupwidth/2 + (2*g-1) * groupwidth / (2*nbars);
    er = errorbar(x, meanDatabySlope(:,g), semDatabySlope(:,g), '.');
    er.Color = [0 0 0];
end %group

xticklabels({'Forward', 'Reverse'})
xlabel('Replay event slope direction')

legend({group(1).name, group(2).name}, 'Location', 'northeastoutside')


end 