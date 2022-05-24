rat = 1;
num = 16;
fs = 2000;

% x = replayEvents(rat).events(num).lfpBuffer(:,1);
% t = ((1:length(x))./2)-400;
% sharpWave = bandpass(x,[5 15], fs);
% ripple = bandpass(x,[130 200],fs);
% gamma = bandpass(x,[20 40],fs);
% 
% figure(1); clf;
% plot(t,[sharpWave,gamma,ripple])
% legend({'Sharp-wave','Gamma','Ripple'})
% xlabel('Time (ms)')
% ylabel('Amplitude (mV)')
% title('Example SWR Components')

figure(2); clf; hold on;
colors = 'rmygbk';
for i = 1:6
    plot(1:4:360,group(1).rat(1).day(1).sleep(1).rip(101).pxn(:,i), colors(i))
end
title('Sample SWR place cell distributions')
xlabel('Positions (in degrees)')
ylabel('Probability')