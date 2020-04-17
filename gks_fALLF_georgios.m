% fALFF analysis and plots

% Go in the folder with the .txt time-courses
cd('/Volumes/Bil04_DataDisk2/Lore/Monica/DATA_APRIL/Preprocessing_new/150_Preprocessed_DREADD/fALFF_handmatig/RsR_R')

% the following will read the names of all .txt files so don't have other
% things there...
d=dir('*.txt');

% some parameters that may need to be adjusted
TR=2;  % TR in seconds
Fs=1/TR; % sampling frequency
L=150;  % number of time points/volumes

% calculate the frequencies
f = Fs*(0:(L/2))/L; 

% this will read the data in matlab
for i=1:numel(d);
    tmp=textread(d(i).name);
    data(:,i)=tmp(:,1);
end

% and this will calculate the FFT
fftdata=fft(data);

% this is to calculate the power
P2=abs(fftdata/L);
P1 = P2(1:L/2+1,:); % the one side of the spectrum (negative frequencies ignored)

% Calculate fALLF
finterest=find(f>=0.05 & f<=0.09);
ftotal=find(f>=0.01 & f<=0.2);
Pinterest=P1(finterest,:);
mPinterest=mean(sqrt(Pinterest),1);
Ptotal=P1(ftotal,:);
mTotal=mean(sqrt(Ptotal),1);
fALLF=mPinterest./mTotal;
mean_fALLF=mean(fALLF);
std_fALLF=std(fALLF);
sem_fALLF=std_fALLF./sqrt(numel(fALLF));

figure; hold on;
%plot all the spectra
plot(f,P1);
% calculate mean and standart error of mean
mSpc=mean(P1,2);
stdSpc=std(P1,1,2);
semSpc=stdSpc./sqrt(size(fftdata,2));
% plot the mean and errorbars
errorbar(f,mSpc,semSpc,'r','LineWidth',3);
title('CNO'); xlabel('frequency (Hz)')

%% now go to the SALINE folder and do the same
cd('/Volumes/Bil04_DataDisk2/Lore/Monica/DATA_APRIL/Preprocessing_new/150_Preprocessed_DREADD/fALFF_handmatig/SALINE/RsR_R');

d=dir('*.txt');
% this will read the data in matlab
for i=1:numel(d);
    tmp=textread(d(i).name);
    data(:,i)=tmp(:,1);
end

% and this will calculate the FFT
fftdata=fft(data);

% this is to calculate the power
P2=abs(fftdata/L);
P1_SAL = P2(1:L/2+1,:); % the one side of the spectrum (negative frequencies ignored)

% Calculate fALLF
finterest=find(f>=0.05 & f<=0.09);
ftotal=find(f>=0.01 & f<=0.2);
Pinterest_SAL=P1_SAL(finterest,:);
mPinterest_SAL=mean(sqrt(Pinterest_SAL),1);
Ptotal_SAL=P1_SAL(ftotal,:);
mTotal_SAL=mean(sqrt(Ptotal_SAL),1);
fALLF_SAL=mPinterest_SAL./mTotal_SAL;
mean_fALLF_SAL=mean(fALLF_SAL);
std_fALLF_SAL=std(fALLF_SAL);
sem_fALLF_SAL=std_fALLF_SAL./sqrt(numel(fALLF_SAL));


figure; hold on;
%plot all the spectra
plot(f,P1_SAL);
% calculate mean and standart error of mean
mSpc_SAL=mean(P1_SAL,2);
stdSpc_SAL=std(P1_SAL,1,2);
semSpc_SAL=stdSpc_SAL./sqrt(size(fftdata,2));
% plot the mean and errorbars
errorbar(f,mSpc_SAL,semSpc_SAL,'k','LineWidth',3);
title('SALINE'); xlabel('frequency (Hz)')

%% Make a ttest and plot the comparison
[~,p]=ttest(fALLF,fALLF_SAL);

figure; hold on;
errorbar(f,mSpc,semSpc,'r','LineWidth',3);
errorbar(f,mSpc_SAL,semSpc_SAL,'k','LineWidth',3);
legend('CNO','SALINE')
title(['fALLF(CNO) = ', num2str(mean_fALLF),' +- ',num2str(sem_fALLF),...
    '/ fALLF(SAL) = ', num2str(mean_fALLF_SAL),' +- ',num2str(sem_fALLF_SAL),...
    ' / p = ',num2str(p)]);


