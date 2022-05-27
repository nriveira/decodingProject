function group = fmr1CircTrack_1_buildDataStruct(dataDir)
% function group = fmr1CircTrack_1_buildDataStruct(dataDir)
%
% PURPOSE:
%  Function to build the main data structure for Emma's circle track FMR1 project.
%
% NOTE:
%  Written for data organization on JBT PC (path below, with rat name as subdir and begin subdirs within that
%
% OUTPUT:
%  group = structure with data for each group (KO vs WT)
%
% JBT
% 11/29/2020
% Colgin Lab

group(1).name = 'WT';
group(2).name = 'KO';

%     *** THIS FUNCTION MUST BE UPDATED FOR EACH NEW RAT ****
group = fmr1CircTrack_0_hardCodeTestData(group); %Add names, dates, theta tet, reward locs
curDir = pwd;

for g = 1:2
    fprintf('%s\n', group(g).name)
    for r = 1:length(group(g).rat)
        fprintf('\tRat %d/%d (%s)\n', r, length(group(g).rat), group(g).rat(r).name);
        for d = 1:length(group(g).rat(r).day)
            fprintf('\t\tDay %d/%d\n', d, length(group(g).rat(r).day));
            cd(strjoin([dataDir '/' group(g).name '/' group(g).rat(r).name '/' group(g).rat(r).day(d).name], ''));
            
            fid = fopen('CellList.txt');
            tmp = textscan(fid, '%s', 'delimiter', '\n');
            unitList = tmp{1};
            fclose(fid);
            
            for b = 1:4
                fprintf('\t\t\tBegin %d\n', b);
                
                group(g).rat(r).day(d).begin(b).dir = strjoin([dataDir '/' group(g).name '/' group(g).rat(r).name '/' group(g).rat(r).day(d).name '/begin' num2str(b)],'');
                cd(group(g).rat(r).day(d).begin(b).dir);
                [post,posy,posx] = LoadCircPos('VT1.nvt');
                
                radPos = circpos(posx,posy); %radial position
                radPos = [post' radPos'];
                
                coords = zeros(length(post), 3);
                coords(:,1) = post;
                coords(:,2) = posx;
                coords(:,3) = posy;
                
                group(g).rat(r).day(d).begin(b).radPos = radPos;
                group(g).rat(r).day(d).begin(b).coords = coords;
                
                tetNums = [];
                uCntr = 0;
                for u = 1:length(unitList)
                    try
                        spkTms = Readtfile(unitList{u});
                        goodRead = 1;
                        uCntr = uCntr + 1;
                    catch
                        goodRead = 0;
                    end
                    
                    if ~isfile(unitList{u})
                        spkTms = [];
                    end
                    
                    if goodRead == 1
                        spkTms = spkTms ./ 10^4;
                        
                        hyphInd = strfind(unitList{u}, '_');
                        tetNum = str2double(unitList{u}(3:hyphInd-1));
                        extInd = strfind(unitList{u}, '.t');
                        clustNum = str2double(unitList{u}(hyphInd+1:extInd-1));
                        
                        group(g).rat(r).day(d).begin(b).unit(uCntr).ID = [tetNum clustNum];
                        group(g).rat(r).day(d).begin(b).unit(uCntr).spkTms = spkTms;
                        
                        tetNums = [tetNums tetNum]; %#ok
                        
                    end
                end
                fprintf('\t\t\t\t%d Units Attached\n', uCntr);
                
                tetNums = unique(tetNums);
                group(g).rat(r).day(d).tetNums = tetNums;
                
                for tt = 1:length(tetNums)
                    cscFn = ['CSC' num2str(tetNums(tt)) '.ncs'];
                    if ~isfile([cscFn(1:end-4) '_LfpForRipDetect'])
                        fprintf('\t\t\t\tAdding LFP for Tetrode #%d\n', tetNums(tt));
                        lfpStruct = read_in_lfp(cscFn);
                        save([cscFn(1:end-4) '_LfpForRipDetect'], 'lfpStruct');
                    end
                        fprintf('\t\t\t\tAdding Wavelet Power for Tetrode #%d\n', tetNums(tt));
                        lfpStruct = read_in_lfp(cscFn);
    
                        wp = get_wavelet_power(lfpStruct.data, lfpStruct.Fs, [1, 250],6);
                        group(g).rat(r).day(d).begin(b).waveletPower(tt).mean = mean(wp,2);
                        group(g).rat(r).day(d).begin(b).waveletPower(tt).std = std(wp,0,2);
                        group(g).rat(r).day(d).begin(b).waveletPower(tt).tetNums = tetNums(tt);
                end
            end %begin            
        end %day
    end
end

cd(curDir); 
end %fnctn
