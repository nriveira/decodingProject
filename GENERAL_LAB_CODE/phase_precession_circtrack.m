function pf = phase_precession_circtrack(spkTms, radPos, pf, lfpStruct)
%
% PURPOSE:
%   Function assesses phase precession, providing spike phases by position and phase-precession stats
%    NOTE: The direction that the rat was running (i.e., clockwise vs. counterclockwise [180 toward 360 vs
%          360 toward 180) will impact the results here. You may have to adjust lines 53-57 (see notes next to those lines)
%
% INPUT:
%     spkTms = vector of spike times for the unit
%     coords = matrix of rat's coordinates by time (output of 'read_in_coords')
%         pf = structure for each place-field (output of 'get_2d_pfs' then passed through 'get_2d_pf_passes')
%  lfpStruct = lfp structure with broad and narrow theta filtered LFPs attached as described in documention
%              for 'get_theta_phase_times'
%
% OUTPUT:
%    pf = inputted structure with addition subfields:
%          1) pf.passSpkTms (spk times for each pass)
%          2) pf.passSpkPhis (spk theta phases for each pass in radians)
%          3) pf.passSpkPos (spk positions in cm, relative to pf start, for each pass)
%          5) pf.passSpkAngPos (angular spk positions, relative to pf start, for each pass)
%          6) pf.passNormSpkPos (within-field spike positions normalized between 0 and 1)
%          7) pf.passPPStats (phiPrecess stats; matrix of slope, intercept, R^2, and pVal for each pass)
%
%
% JB Trimper
% 09/2019
% Colgin Lab



angToDegScale = (pi*100)/360; %because track has 100cm diameter;


for p = 1:length(pf)
    pf(p).passSpkTms = cell(size(pf(p).passes,1),1); %cell containing spkTms for each pass
    pf(p).passSpkPhis = cell(size(pf(p).passes,1),1); %cell containing theta spk phase for each pass
    pf(p).passSpkPos = cell(size(pf(p).passes,1),1); %cell containing spkTms for each pass
    pf(p).passSpkAngPos = cell(size(pf(p).passes,1),1); %cell containing angular spike position relative to field start
    pf(p).passNormSpkPos = cell(size(pf(p).passes,1),1); %cell containing within-field spike positions normalized between 0 and 1
    pf(p).passPPStats = nan(size(pf(p).passes,1),4); %matrix of slope, intercept, R^2, and pVal for each pass
end


phiTms = get_theta_phase_times(lfpStruct); %get theta phase times (pk, fa, tr, ri)
phiVctr = get_asym_theta_phi_vector(lfpStruct); %get a theta phase vector, interpolated from phase times


for p = 1:length(pf) %for each place-field
    if ~isempty(pf(p).passes) %if there were passes through that field
        
        
        for pp = 1:size(pf(p).passes,1) %for each pass
            
            
            %If rat was running from 360 toward 180, use this line
            startPos = pf(p).radPos(end);
            
            %If rat was running from 180 toward 360, use this line
            %             startPos = pf(p).radPos(1);
            
            
            startInd = pf(p).passes(pp,1);
            startTm = radPos(startInd,1); %start time for the pass
            
            endInd = pf(p).passes(pp,2);
            endTm = radPos(endInd,1); %end time for the pass
            
            passSpks = spkTms(spkTms>=startTm & spkTms<=endTm); %spikes during the pass
            
            if length(passSpks)>=5 %If at least 5 spikes were present
                
                firstCycle = find(phiTms(1,:)<=passSpks(1), 1, 'Last'); %first theta cycle of pass
                lastCycle = find(phiTms(1,:)>=passSpks(end), 1, 'First'); %last theta cycle of pass
                
                if range(firstCycle:lastCycle)+1 >= 4 %and those spikes covered a range of at least 4 theta cycles
                    
                    spkPhis = zeros(1,length(passSpks)); %pre-allocated for spike phases
                    absSpkPos = zeros(1,length(passSpks)); %pre-allocated for spike's angular position relative to whole track
                    
                    for st = 1:length(passSpks) %for each spike
                        
                        spkLfpInd = find(lfpStruct.ts<=passSpks(st), 1, 'Last'); %get the closest LFP ind
                        spkPhis(st) = phiVctr(spkLfpInd); %get the theta phase at which spike occurred
                        
                        posInd = find(radPos(:,1)<=passSpks(st), 1, 'Last'); %closest position index
                        absSpkPos(st) = radPos(posInd,2); %radial position at which spike occurred
                        
                    end %spkTms
                    
                    spkAngPos = abs(rad2deg(circ_dist(deg2rad(absSpkPos), deg2rad(startPos)))); %distance in degrees from spike position to field start position
                    spkPos = spkAngPos .* angToDegScale; %convert angular distances of spike positions from peak to centimeter distances
                    
                    % Get position of spikes normalized from 0 to 1
                    numSpks = st;
                    radFldPos = pf(p).radPos';
                    tmp = rescale([absSpkPos'; radFldPos], 0, 1);
                    normRadSpkPos = tmp(1:numSpks)';
                    normRadSpkPos = 1 - normRadSpkPos; %because the rats ran 360 to 0 instead of 0 to 360

                    [beta,R2,pVal] = CircularRegression(spkPos,spkPhis); %calculate circular/linear regression of spike phases/positions
                    
                    %store information in output structure
                    pf(p).passSpkTms{pp} = passSpks';
                    pf(p).passSpkPhis{pp} = spkPhis;
                    pf(p).passSpkAngPos{pp} = spkAngPos;
                    pf(p).passSpkPos{pp} = spkPos;
                    pf(p).passNormSpkPos{pp} = normRadSpkPos; 
                    pf(p).passPPStats(pp,:) = [beta R2 pVal]; %slope intercept R^2 sigVal
                    
                end %if enough theta cycles covered
            end %if enough spks
            
            
        end %pass
    end %if passes
end %placefields






end %fnctn