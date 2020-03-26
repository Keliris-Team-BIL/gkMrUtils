function gk_alias_range(BPM_FromTo,TRsec,critF)
% Usage: gk_alias_range(BPM_FromTo,TRsec,critF)
%
% Function to calculate the alias frequency for a range of signal
% frequencies and identify and plot those that would give be below a 
% threshold frequency e.g. 0.1 or 0.2 for resting state analysis
%
% INPUT
%   BPM_FromTo  : a vector with two elements defining the freq range in BPM
%   TRsec       : the sampling interval (TR) in seconds
%   critF       : the critical frequency in Hz
%
% EXAMPLE
%   gk_alias_range([0 400],0.6,0.2)
%  
% GAK, May 2019



rangeBPM=BPM_FromTo(1):0.01:BPM_FromTo(2);
rangeHz=rangeBPM/60.0;
for i=1:numel(rangeHz)
    aliased_freq(i)=min(abs(rangeHz(i)-[1:10]./TRsec));
end

goodFreqInd=find(aliased_freq>critF);
freqs=zeros(size(aliased_freq));
freqs(goodFreqInd)=0.5/TRsec;
badFreqStart=find(diff(goodFreqInd)>1);

figure; hold on;
area(rangeHz*60,freqs,'Facecolor',[0.7 1 0.7],'LineStyle','none');
plot(rangeHz*60,aliased_freq,'r'); hold on;
plot([0 rangeHz(end)*60],[critF critF],'k--')    
ylabel('Aliased Frequency')
xlabel('Signal -> BPM')
xlim([rangeHz(1) rangeHz(end)]*60);
ylim([0 0.5/TRsec])

fprintf('Avoid frequency ranges: \n');
try
    fprintf('\t%.1f - %.1f Hz\n',60*rangeHz([1 goodFreqInd(1)-1]));
end
for i=1:numel(badFreqStart)
    fprintf('\t%.1f - %.1f Hz\n',60*rangeHz(goodFreqInd(badFreqStart(i):badFreqStart(i)+1)));
end
