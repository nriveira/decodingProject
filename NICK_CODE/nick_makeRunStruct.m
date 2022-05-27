function runEvents = nick_makeRunStruct(group, saveDir)
    %saveDir = "C:/Users/nick/Projects/DATA_STRUCTS";
    for g = 1:length(group)
        fprintf('Group %d\n', g);
        for r = 1:length(group(g).rat)
            fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);
            for d = 1:length(group(g).rat(r).day)
                fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day));
                runEvents(g).name = group(g).name;
                for b = 1:length(group(g).rat(r).day(d).begin)
                    curDir = pwd;
                    beginDir = group(g).rat(r).day(d).begin(b).dir;
                    cd(beginDir)
                    lfpDirs = dir('CSC*LfpForRipDetect.mat');
    
                    load(lfpDirs(1).name);
                    dataCombined = zeros(length(lfpStruct.data), length(lfpDirs));
                    dataTs = lfpStruct.ts;
                    Fs = lfpStruct.Fs;
                    sigThetaInd = 1;
    
                    for lfp = 1:length(lfpDirs)
                        load(lfpDirs(lfp).name);
                        dataCombined(:,lfp) = lfpStruct.data;
                    end
                    normFactorDayMean = group(g).rat(r).day(d).dayBeginNormMean;
                    normFactorDaySTD = group(g).rat(r).day(d).dayBeginNormSTD;
                    normFactorBeginMean = group(g).rat(r).day(d).begin(b).beginNorm;
                    normFactorBeginSTD = group(g).rat(r).day(d).begin(b).beginSTD;
                    
                    [~, edgeInds, ~] = find_run_bouts(group(g).rat(r).day(d).begin(b).coords,5,5);
                    for ind = 1:length(edgeInds)
                        signal = dataCombined(edgeInds(ind,1):edgeInds(ind,2),:);
                        [~, cutSignal] = cut2theta(signal, Fs);
                        if(~isempty(cutSignal))
                            runEvents(g).events(sigThetaInd).cutTheta = cutSignal;
    
                            normWaveletPowerDay = zeros(250, length(cutSignal), size(cutSignal,2));
                            normWaveletPowerBegin = zeros(250, length(cutSignal), size(cutSignal,2));
            
                            % Manually z-score
                            for tet = 1:size(cutSignal,2)
                                rawWaveletPower = get_wavelet_power(cutSignal(:,tet),2000,[1,250],6);
                                for t = 1:length(cutSignal)
                                    normWaveletPowerDay(:,t,tet) = (rawWaveletPower(:,t)-normFactorDayMean(tet,:)')./normFactorDaySTD(tet,:)';
                                    normWaveletPowerBegin(:,t,tet) = (rawWaveletPower(:,t)-normFactorBeginMean(tet,:)')./normFactorBeginSTD(tet,:)';
                                end
                            end
                            runEvents(g).events(sigThetaInd).wavelet_power_day = mean(normWaveletPowerDay,3);
                            runEvents(g).events(sigThetaInd).wavelet_power_begin = mean(normWaveletPowerBegin,3);
                            sigThetaInd = sigThetaInd+1;
                        end
                    end                   
                end %begin
            end %day
        end %rat
    end %group

    cd(saveDir); save('runDataStruct','runEvents','-v7.3')
end %function