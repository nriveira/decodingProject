function quickNormRm_circTrack(varagin)
% function quickNormRm_circTrack(dayFolder)
%
% PURPOSE:
%    To plot ratemaps for each unit from one day to check if cells are
%    spatially modulated.
%
% INPUT:
%    dayFolder (optional) = directory address for folder from this day
%           ex. dayFolder = 'E:\FMR1_SWRS\RAW_DATA\Rat326Z\2021-06-11';
%           If not included, function uses current working directory.
%
% OUTPUT:
%    Figures. Also saved as .pngs in day folder.
%
%  NOTE:
%    Folder must contain "TTList.txt" with accepted clusters.
%
%    All ratemaps are normalized to max firing rate in that day. See
%    function quickRawRm for ratemaps normalized within each begin window.
%
% MMD
% 06/2021
% Colgin Lab

%% INITIALIZE

if nargin == 0
    dayFolder = pwd;
end %use input or not

spatBinSz = 4; %degrees
velFilt = 1;
durCrit = 1;

gWinDur = 30; %degrees
gWinStd = gWinDur/2; %degrees

numBins = 360/spatBinSz; %360/5 degree bins

dashInd = strfind(dayFolder, '\');
dashInd = dashInd(end);
dayName = dayFolder(dashInd+1:end);

colMap = [.5 .5 .5; jet(length(-0.01:0.01:1))]; %jet with grey nan option

cd(dayFolder)
clustlist = ReadFileList('TTList.txt');

% spOpts = [1:4; 5:8; 9:12; 13:16; 17:20; 21:24]; %subplot options
spOpts = [1:5; 6:10; 11:15; 16:20; 21:25; 26:30]; %subplot options


allRms = nan(numBins, 4, length(clustlist));
uMaxFr = nan(1,length(clustlist)); %for adding to plot

%% GET DATA

for u = 1:length(clustlist)
    fprintf('\tUnit %d/%d\n', u, length(clustlist))
    unitRms = nan(numBins, 4);
    for b = 1:4
        
        cd(['begin' num2str(b)])
        
        filename = clustlist{u};
        if ~isfile(filename)
            fprintf('\t\t%s - Begin %d\n', filename, b)
            error('Need to write files for this tetrode!') %I forget to do this a lot
        end
        spkTms = Readtfile(filename);
        spkTms = spkTms ./ 10^4; %convert to seconds
        
        [radPos, coords] = read_in_circtrack_coords('VT1.nvt');
        
        rateMap = get_ratemap_circtrack(spkTms, coords, radPos, spatBinSz, velFilt, durCrit);
        smRateMap = smooth_circtrack_ratemap(rateMap, gWinDur, gWinStd);
        
        
        unitRms(:,b) = smRateMap;
        cd ../
    end %begin
    
    maxFr = max(unitRms(:));
    uMaxFr(u) = maxFr;
    normUnitRms = unitRms ./ maxFr;
    
    allRms(:,:,u) = normUnitRms;
    
end %unit

%% MAKE FIGS

numFigs = ceil(length(clustlist)/6); %how many figure do we need - 6 units per fig

uCntr = 1;
for f = 1:numFigs
    figTitle = [dayName '_NormRateMaps_Figure' num2str(f)];
    
    figure ('Position', [680 63 832 915], 'Name', figTitle);
    
    for figUCntr = 1:6
        if uCntr <= length(clustlist)
            uName = clustlist{uCntr};
            undInd = strfind(uName, '_');
            uName = [uName(1:undInd-1) '\' uName(undInd:end)];
            
            subplot(6,5, spOpts(figUCntr,1))
            text(0 ,0.5, {uName; ['Max = ' num2str(round(uMaxFr(uCntr),1)) ' Hz']});
            axis off
            
            for b = 1:4
                subplot(6,5, spOpts(figUCntr,b+1))
                
                plotMap = allRms(:,b,uCntr);
                CirHeatmap({plotMap}, 'CircType', 'o');
                
                colormap(colMap)
                caxis([-0.01 1])
                
                if figUCntr == 6 || uCntr == length(clustlist) %if last unit in this fig
                    xlabel('Position (cm)');
                end
                
                if figUCntr == 1
                    title(['Begin ' num2str(b)]);
                end %if first unit
                
                
            end %begin
            uCntr = uCntr + 1;
        end %if there are still units to plot
    end %units in this fig
    
    saveas(gcf, figTitle, 'png')  %save png of this fig in day folder
    
end %figures


end %function