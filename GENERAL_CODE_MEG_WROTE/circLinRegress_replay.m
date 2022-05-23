function [xAx, com, calphase, slope, r2, p] = circLinRegress_replay(pxn, radBinCtrs, timeAx, bins2use)
% function [xAx, com, calphase, slope, r2, p] = circLinRegress_replay(pxn, radBinCtrs, timeAx, bins2use)
%
% NOTE:
%   See Cir_reg function from CZ for reference. This function is edited
%   slightly from original for clarity and comments. Additional changes
%   noted in function.
%
% PURPOSE:
%   Perform the circular regression for replay events on the circle track.
%
% INPUT:
%   pxn = decoded probability distribution across replay event.
%   radBinCtrs = angle in radians of the bin centers for the pxn.
%   timeAx = time in seconds of the time bins for the pxn.
%   bins2use = which time bins from pxn to use for regression; bins where
%       there were no spikes (aka pxn is chance across all position bins)
%       should be excluded.
%
% OUTPUT:
%   com = center of mass of the distribution for each time bin (with bins
%       removed where there were no spikes)
%   r2 = r^2 value for the regression, drawn from 1000 surrogate r2 values
%   calphase = calculated values for the regression line across the time
%       axis (with bins removed where there were no spikes)
%   xAx = x values for the regression (timeAx with bins removed that were
%       not used because there were no spikes)
%   p = p-value for the regression significance
%   slope = slope of the regression line, in rad/s
%
% MMD
% Colgin Lab
% 03/2022

%% CHECK

if size(radBinCtrs,1) == 1
    radBinCtrs = radBinCtrs';
end %align correctly

if length(radBinCtrs) ~= size(pxn,1)
    error('Position bins must be the same for pxn and bin centers')
end

%% INITIALIZE

% bound = round(1/timeAx(end)); %max of entire track apprx. (orignial code = 8)

drawnum = 1000; %how many surrogate r2 values to get

%% CALCULATE

com = NaN(1,size(pxn,2)); % initialize center of mass
for t = 1:length(bins2use) %across time bins to use
    bin = bins2use(t); %which time bin it is
    
    com(bin) = circ_mean(radBinCtrs, pxn(:,bin)); %weighted circular mean based on prob distribution
end %bins - i

com = wrapTo2Pi(com); %now we have the center of mass of each time bin

com = com(bins2use);
xAx = timeAx(bins2use);

% [r2, p] = circ_corrcl(com, xAx);
%
% if r2 >= 0.5

[beta, ~, p] = CircularRegression(xAx, com);
slope = beta(1);

calphase = beta(1)*xAx + beta(2);
% calphase = wrapTo2Pi(calphase);

% else
%    slope = NaN;
%    calphase = NaN;
% end

%
% if r2 < 0
%     keyboard
% end
% 
% para = circ_lin_regress(xAx, com, bound);
% 
% calphase = 2*pi*para(1,1)*xAx + para(1,2);
% tmpSlope = para(1,1)*2*pi;
% 
% if round(tmpSlope,3) ~= round(slope,3)
%     keyboard
% end

% tmpCheck
%
% xAx = xAx';
% com = com';
% 
% yRes = zeros(drawnum,length(xAx));
% yTot = zeros(drawnum,length(xAx));
% for i = 1:length(xAx)
%     yFit = repmat(calphase(i), drawnum, 1);
%     yRes(:,i) = circdistance(com(i), yFit, 1);
%     yTot(:,i) = circdistance(com(i), circ_mean(com), 1);
% end %time bins
% 
% SSres = sum(yRes(:));
% SStot = sum(yTot(:));
% r2 = 1 - SSres/SStot;
%
% yRes = zeros(drawnum,length(xAx));
% yTot = zeros(drawnum,length(xAx));
% for i = 1:length(xAx)
%     yFit = calphase(i);
%     yRes(:,i) = circ_dist(com(i), yFit)^2;
%     yTot(:,i) = circ_dist(com(i), circ_mean(com'))^2;
% end %time bins
% 
% SSres = sum(yRes(:));
% SStot = sum(yTot(:));
% r2 = 1 - SSres/SStot;
% if r2<0
%     r2 = NaN;
% end

r2 = circLinRegress_r2(com, calphase);

end %function