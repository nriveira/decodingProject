function same_caxis(useMin, varagin)
% function same_caxis(useMin, figHandle)
%  
% PURPOSE:
%   To make the colobar axis limits identical in all subplots in a figure.
% 
% INPUT:
%   useMin = use the minimum on the figure for min caxis.
%   figHandle (optional) = handle for figure to edit. If not included, code
%       will edit most recently created figure.
% 
% MMD
% 12/2021
% Colgin Lab

if nargin < 2
    figHandle = gcf;
end %use input or not
% 
% if nargin == 1
%     useMin = 0;
% end %use input or not

% numSub = numel(figHandle.Children); %get number of subplots
figure(figHandle)
findAx = findall(gcf, 'type', 'axes');
numSub = numel(findAx);
for n = 1:numSub
    pos1(n) = round(findAx(n).Position(1),2); %rows
    pos2(n) = round(findAx(n).Position(2),2); %columns
end %all subplots
nCols = numel(unique(pos1));
nRows = numel(unique(pos2));

cMin = inf; %initialize
cMax = -inf;

for sp = 1:numSub
    subplot(nRows, nCols, sp)
    ax = gca;
    
    if cMax < ax.CLim(2)
        cMax = ax.CLim(2);
    end %max
    if cMin > ax.CLim(1)
        cMin = ax.CLim(1);
    end %min
end %subplots to check  min and max

if useMin == 0
    cMin = 0;
end %set min as 0

for sp = 1:numSub
    subplot(nRows, nCols, sp)
    caxis([cMin cMax])
end %subplot to change axes

end %function