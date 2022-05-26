function group = fmr1CircTrack_7_normalizeWP(group)
    for g = 1:length(group)
        fprintf('Group %d\n', g);
        for r = 1:length(group(g).rat)
            fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);
            for d = 1:length(group(g).rat(r).day)
                fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day));
                
                normMean = zeros(length(group(g).rat(r).day(d).tetNums), 250, length(group(g).rat(r).day(d).begin));
                normSTD = zeros(length(group(g).rat(r).day(d).tetNums), 250, length(group(g).rat(r).day(d).begin));
                for b = 1:length(group(g).rat(r).day(d).begin)
                    fprintf('\t\t\tBegin %d/%d\n', b, length(group(g).rat(r).day(d).begin));
                    for tt = 1:length(group(g).rat(r).day(d).tetNums)
                        normMean(tt,:,b) = group(g).rat(r).day(d).begin(b).waveletPower(tt).mean;
                        normSTD(tt,:,b) = group(g).rat(r).day(d).begin(b).waveletPower(tt).std.^2;
                    end
                    group(g).rat(r).day(d).begin(b).beginNorm = normMean(:,:,b);
                    group(g).rat(r).day(d).begin(b).beginSTD = sqrt(normSTD(:,:,b));
                end
                group(g).rat(r).day(d).dayBeginNormMean = mean(normMean,3);
                group(g).rat(r).day(d).dayBeginNormSTD = sqrt(mean(normSTD,3));

                normMean = zeros(length(group(g).rat(r).day(d).tetNums), 250, length(group(g).rat(r).day(d).sleep));
                normSTD = zeros(length(group(g).rat(r).day(d).tetNums), 250, length(group(g).rat(r).day(d).sleep));
                for s = 1:length(group(g).rat(r).day(d).sleep)
                    fprintf('\t\t\tSleep %d/%d\n', s, length(group(g).rat(r).day(d).sleep));
                    for tt = 1:length(group(g).rat(r).day(d).tetNums)
                        normMean(tt,:,s) = group(g).rat(r).day(d).sleep(s).waveletPower(tt).mean;
                        normSTD(tt,:,s) = group(g).rat(r).day(d).sleep(s).waveletPower(tt).std.^2;
                    end
                    group(g).rat(r).day(d).sleep(s).sleepNorm = normMean(:,:,s);
                    group(g).rat(r).day(d).sleep(s).sleepSTD = sqrt(normSTD(:,:,s));
                end
                group(g).rat(r).day(d).daySleepNormMean = mean(normMean,3);
                group(g).rat(r).day(d).daySleepNormSTD = sqrt(mean(normSTD,3));

            end
        end
    end
end