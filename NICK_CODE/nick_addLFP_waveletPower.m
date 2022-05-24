for g = 1:length(group)
    fprintf('Group %d\n', g);
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day));

            sri = 1; % Significant replay index
            rep_max = 1;
            rep_min = 100000;
            replayEvents(g).name = group(g).name;
            
            for s = 1:length(group(g).rat(r).day(d).sleep)
                curDir = pwd; 
                sleepDir = group(g).rat(r).day(d).sleep(s).dir;
                cd(sleepDir)
                lfpDirs = dir('*.mat');
                lfpStructs = [];

                load(lfpDirs(1).name);
                dataCombined = zeros(length(lfpStruct.data), length(lfpDirs));
                dataTs = lfpStruct.ts;

                for lfp = 1:length(lfpDirs)
                    load(lfpDirs(lfp).name);
                    lfpStructs = [lfpStructs; lfpStruct];
                    dataCombined(:,lfp) = lfpStruct.data;
                end

                for u = 1:length(group(g).rat(r).day(d).sleep(s).popEv)
                    if(group(g).rat(r).day(d).sleep(s).popEv(u).r2 >= 0.5)
                        times = group(g).rat(r).day(d).sleep(s).popEv(u).tms;
                        beginInd = find(dataTs>=times(1),1);
                        endInd = find(dataTs>=times(2),1);
                        swrLength = endInd-beginInd;

                        if(swrLength < rep_min)
                            rep_min = swrLength;
                        elseif(swrLength > rep_max)
                            rep_max = swrLength;
                        end

                        % Pre/post event indices (using 0.4 * 2000 = 800
                        % points)
                        preBeginInd = max(1,beginInd-800);
                        postEndInd = min(length(dataCombined), endInd+800);

                        replayEvents(g).events(sri).slope = group(g).rat(r).day(d).sleep(s).popEv(u).slope;
                        replayEvents(g).events(sri).r2 = group(g).rat(r).day(d).sleep(s).popEv(u).r2;
        
                        replayEvents(g).events(sri).lfpSWR = dataCombined(beginInd:endInd,:);
                        [replayEvents(g).events(sri).welch, replayEvents(g).freq] = pwelch(replayEvents(g).events(sri).lfpSWR, [],[],[], 2000);
                        replayEvents(g).events(sri).wavelet_power = get_wavelet_power(dataCombined(beginInd:endInd,:),2000,[25,250],6,0,1);
    
                        replayEvents(g).events(sri).lfpBuffer = dataCombined(preBeginInd:postEndInd,:);
                        replayEvents(g).events(sri).wavelet_powerPre = get_wavelet_power(replayEvents(g).events(sri).lfpBuffer,2000,[25,250],6,0,1);
                        sri = sri+1;
                    end
                end
                cd(curDir)
            end
            replayEvents(g).max_length = rep_max;
            replayEvents(g).min_length = rep_min;
        end
    end
end