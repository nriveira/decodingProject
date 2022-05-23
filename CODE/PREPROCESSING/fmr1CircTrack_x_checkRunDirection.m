function fmr1CircTrack_x_checkRunDirection(group)
% function fmr1CircTrack_x_checkRunDirection(group)
% 
% PURPOSE:
%   For visual inspection of rad pos for the first begin of first day to 
%   check run direction.
% 
% INPUT:
%   group struct
% 
% OUTPUT:
%   Figures, not saved.
% 
% MMD
% 6/2021
% Colgin Lab


d = 1; %assuming same direction all days
b = 1;

for g = 1:2
    for r = 1:length(group(g).rat)
    radPos = group(g).rat(r).day(d).begin(b).radPos;
    
    figure
    plot(radPos((1:100:end),2))
    title([group(g).rat(r).name])
    ylim([0 360])    
    end %rat
    
    
end %group


end %function