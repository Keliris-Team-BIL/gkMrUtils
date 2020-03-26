function gk_alias(signalBPM,TRsec,scanDuration)
% Usage: gk_alias(signalBPM,TRsec,scanDuration)
%
% Function to calculate the expected alias frequency for a signal
% and plot it for a specific scan duration (shorter => broader peaks) 
%
% INPUT
%   signalBPM   : the frequency of the signal (e.g. respiration) in BPM
%   TRsec       : the sampling interval (TR) in seconds
%   scanDuration: the scan duration in seconds
%
% EXAMPLE
%   gk_alias(70,0.6,300)
%  
% GAK, May 2019

fresp_bpm=signalBPM; %bpm
fresp=fresp_bpm/60.0;

fs=5000;
dt=1/fs;
T=scanDuration; % signal duration in seconds
t=0:dt:T;

TR=TRsec;
vol=rand/100:TR:T;

Sresp=sin(2*pi*fresp*t);

Sresp_alias=sin(2*pi*fresp*vol);

% if scan duration longer than 30 sec only plot the first 30 sec.
plot_dur=min(scanDuration,60); 
figure;
subplot(2,1,1); hold on;
plot(t(t<=plot_dur),Sresp(t<=plot_dur));
ylabel('Signal')
xlabel('Time -> seconds')
plot(vol(vol<=plot_dur),Sresp_alias(vol<=plot_dur),'ro');
ylim([-1.1 1.1]);
subplot(2,1,2);
spc=fft(Sresp_alias);
P2=2*abs(spc/numel(spc));
P1=P2(1:floor(numel(spc)/2)+1);
f=linspace(0,0.5/TR,numel(P1));
plot(f,P1,'r');
ylabel('Power')
xlabel('Frequency -> Hz')
aliased_freq=abs(fresp-[1:10]./TR);
fprintf("Expected: %.3f\n",min(aliased_freq))