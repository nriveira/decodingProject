function replayEvents = nick_makeReplayStruct(group, saveDir)
    %saveDir = "C:/Users/nick/Projects/DATA_STRUCTS";
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

                tic
                for s = 1:length(group(g).rat(r).day(d).sleep)
                    curDir = pwd; 
                    sleepDir = group(g).rat(r).day(d).sleep(s).dir;
                    cd(sleepDir)
                    lfpDirs = dir('CSC*LfpForRipDetect.mat');
    
                    load(lfpDirs(1).name);
                    dataCombined = zeros(length(lfpStruct.data), length(lfpDirs));
                    dataTs = lfpStruct.ts;
    
                    for lfp = 1:length(lfpDirs)
                        load(lfpDirs(lfp).name);
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
                            signal = dataCombined(beginInd:endInd,:);
                            normWaveletPowerDay = zeros(250, length(signal), size(signal,2));
                            normWaveletPowerSleep = zeros(250, length(signal), size(signal,2));
        
                            normFactorDayMean = group(g).rat(r).day(d).daySleepNormMean;
                            normFactorDaySTD = group(g).rat(r).day(d).daySleepNormSTD;
                            normFactorSleepMean = group(g).rat(r).day(d).sleep(s).sleepNorm;
                            normFactorSleepSTD = group(g).rat(r).day(d).sleep(s).sleepSTD;
    
                            % Manually z-score
                            for tet = 1:size(signal,2)
                                rawWaveletPower = get_wavelet_power(signal(:,tet),2000,[1,250],6);
                                for t = 1:length(signal)
                                    normWaveletPowerDay(:,t,tet) = (rawWaveletPower(:,t)-normFactorDayMean(tet,:)')./normFactorDaySTD(tet,:)';
                                    normWaveletPowerSleep(:,t,tet) = (rawWaveletPower(:,t)-normFactorSleepMean(tet,:)')./normFactorSleepSTD(tet,:)';
                                end
                            end
    
                            replayEvents(g).events(sri).wavelet_power_day = mean(normWaveletPowerDay,3);
                            replayEvents(g).events(sri).wavelet_power_sleep = mean(normWaveletPowerSleep,3);
     
                            normWaveletPowerDayPre = zeros(250, length(signal), size(signal,2));
                            normWaveletPowerSleepPre = zeros(250, length(signal), size(signal,2));
    
                            replayEvents(g).events(sri).lfpBuffer = dataCombined(preBeginInd:postEndInd,:);
                            signal = dataCombined(preBeginInd:postEndInd,:);
    
                             % Manually z-score
                            for tet = 1:size(signal,2)
                                rawWaveletPowerPre = get_wavelet_power(signal(:,tet),2000,[1,250],6);
                                for t = 1:length(signal)
                                    normWaveletPowerDayPre(:,t,tet) = (rawWaveletPowerPre(:,t)-normFactorDayMean(tet,:)')./normFactorDaySTD(tet,:)';
                                    normWaveletPowerSleepPre(:,t,tet) = (rawWaveletPowerPre(:,t)-normFactorSleepMean(tet,:)')./normFactorSleepSTD(tet,:)';
                                end
                            end
                            replayEvents(g).events(sri).wavelet_power_dayPre = mean(normWaveletPowerDayPre,3);
                            replayEvents(g).events(sri).wavelet_power_sleepPre = mean(normWaveletPowerSleepPre,3);
    
                            sri = sri+1;
                        end
                    end
                    toc
                    cd(curDir)
                end % sleep
                replayEvents(g).max_length = rep_max;
                replayEvents(g).min_length = rep_min;
            end % day
        end % rat
    end % group

    cd(saveDir); save('replayDataStruct','replayEvents','-v7.3')
end % function