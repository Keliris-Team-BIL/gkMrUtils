%% This script contains the analysis of the rsfMRI scans after normalization in SPM => KORD study
% Data structure can be found in file 'Analysis_KORD_datastructure'

% subtract fitted (linear + exponential fit) sham timecourse from KORD timecourses per voxel

rawdatapath='/Users/gkeliris/Documents/DATA/Lore_KORD';
% Define the data structure
dat=KORD_datapaths(rawdatapath);

% Define output path for the analysis
processedDataPath=cat(2,rawdatapath,'/ProcessedData_zscore');
if ~isfolder(processedDataPath)
    mkdir(processedDataPath);
    mkdir(processedDataPath,'/SHAM');
    mkdir(processedDataPath,'/KORD');
    mkdir([processedDataPath,'/SHAM'],'LOW_CONC');
    mkdir([processedDataPath,'/KORD'],'LOW_CONC');
    mkdir([processedDataPath,'/KORD'],'HIGH_CONC');
end

% Get voxel dimensions, timepoints from 1st sham scan
scaninfo=niftiinfo(dat.sham.lowconc.raw{1});
N_timepoints    =scaninfo.ImageSize(4);
N_slices        =scaninfo.ImageSize(3);
Vox_x           =scaninfo.PixelDimensions(1);
Vox_y           =scaninfo.PixelDimensions(2);
Vox_z           =scaninfo.PixelDimensions(3);

% Define the analysis parameters
baseline=1:20;
% spcSmooth=2*[Vox_x Vox_y Vox_x];   % Use Vox_x for the slices as well to make smoothing isometric
spcSmooth=3*Vox_x;   % Use Vox_x for the slices as well to make smoothing isometric


%% PRE-PROCESSING (can be skipped if already done)
%--------------------------------------------------------------------------
to_pch=1; to_zscr=2;
normtype=to_zscr;
% SHAM LOW-CONC
gk_kord_preprocess(dat.sham.lowconc.normalized,[processedDataPath,'/SHAM/LOW_CONC/'],...
    dat.ROIs.mask{1},spcSmooth,baseline,normtype);
% KORD LOW-CONC
gk_kord_preprocess(dat.kord.lowconc.normalized,[processedDataPath,'/KORD/LOW_CONC/'],...
    dat.ROIs.mask{1},spcSmooth,baseline,normtype);
% KORD HIGH-CONC
gk_kord_preprocess(dat.kord.highconc.normalized,[processedDataPath,'/KORD/HIGH_CONC/'],...
    dat.ROIs.mask{1},spcSmooth,baseline,normtype);


%% AVERAGE SHAMs TO GET A SINGLE TIME SERIES (can be skipped if already done)
%--------------------------------------------------------------------------
to_pch=1; to_zscr=2;
normtype=to_zscr;
if normtype==1
    shamFiles=dir([processedDataPath,'/SHAM/LOW_CONC/pch_mSmth_twu*.nii']);
elseif normtype==2
    shamFiles=dir([processedDataPath,'/SHAM/LOW_CONC/zscr_mSmth_twu*.nii']);
end
shams_to_average=[1 2 4 6];
info=niftiinfo(fullfile(shamFiles(1).folder,shamFiles(1).name));
for i=1:numel(shams_to_average)
    n=shams_to_average(i);
    m(:,:,:,:,i)=niftiread(fullfile(shamFiles(n).folder,shamFiles(n).name));
end
mean_sham=mean(m,5);
niftiwrite(mean_sham,fullfile(shamFiles(1).folder,'mean_sham.nii'),info);

%% Fit linear curve (volume 5:20) and double exponential curve (volume 25:120) to average sham scan (per voxel)
%--------------------------------------------------------------------------
ft = gk_timecourseFit(fullfile(shamFiles(1).folder,'mean_sham.nii'),dat.ROIs.mask{1},false);

%% subtract fitted curves from KORD scans
%--------------------------------------------------------------------------
%% Calculate the AUC values to perform statistics on
aucMapKL = gk_calc_AUC(dat,[processedDataPath,'/SHAM/LOW_CONC/'],...
    [processedDataPath,'/KORD/LOW_CONC/'],0.95,[80:100],normtype);
aucMapKH = gk_calc_AUC(dat,[processedDataPath,'/SHAM/LOW_CONC/']...
    ,[processedDataPath,'/KORD/HIGH_CONC/'],0.95,[80:100],normtype);



%% Calculate the laterality using AveLI
%------------------------------------------------------
%% OPTIONAL TO START FROM THIS PART
rawdatapath='/Users/gkeliris/Documents/DATA/Lore_KORD';
processedDataPath=cat(2,rawdatapath,'/ProcessedData_zscore');

% Define the data structure
dat=KORD_datapaths(rawdatapath,0);


% For positive ROIs
kordFiles=dir([processedDataPath,'/KORD/HIGH_CONC/auc_z_adj_zscr_mSmth_twu*.nii']);
kordFilesCell=cellfun(@(x,y) fullfile(x,y), {kordFiles(:).folder},{kordFiles(:).name}, 'Uniformoutput',false);

res = gk_aveLI(kordFilesCell, dat.ROIs.left([2:6 9:11]), dat.ROIs.right([2:6 11:13]), 1);

% For negative ROIs (ACA, TH, VC)
kordFiles=dir([processedDataPath,'/KORD/HIGH_CONC/neg_auc_z_adj_zscr_mSmth_twu*.nii']);
if isempty(kordFiles)
  kordFiles=dir([processedDataPath,'/KORD/HIGH_CONC/auc_z_adj_zscr_mSmth_twu*.nii']);
  for i=1:numel(kordFiles)
      data_info   = niftiinfo(fullfile(kordFiles(i).folder,kordFiles(i).name));
      data        = niftiread(data_info);
      data2       = -1*data;
      niftiwrite(data2,fullfile(kordFiles(i).folder,['neg_',kordFiles(i).name]),data_info);
  end
  kordFiles=dir([processedDataPath,'/KORD/HIGH_CONC/neg_auc_z_adj_zscr_mSmth_twu*.nii']);
end
kordFilesCell=cellfun(@(x,y) fullfile(x,y), {kordFiles(:).folder},{kordFiles(:).name}, 'Uniformoutput',false);

res = gk_aveLI(kordFilesCell, dat.ROIs.left([1 7:8]), dat.ROIs.right([1 7:8]), 1);





%% THE THINGS BELOW ARE FOR PLOTTING THE TIME COURSES
%-----------------------------------------------------
%% Header: Define file paths
rawdatapath='/Users/gkeliris/Documents/DATA/Lore_KORD';
processedDataPath=cat(2,rawdatapath,'/ProcessedData_zscore');
statsFile='/Users/gkeliris/Documents/DATA/Lore_KORD/Zscore_smth3_95_80-100/Highconc_spmT_0001.nii';
zScoreTcFilesHigh=dir([processedDataPath,'/KORD/HIGH_CONC/z_adj_zscr_mSmth_twu*.nii']);
zScoreTcFilesLow=dir([processedDataPath,'/KORD/LOW_CONC/z_adj_zscr_mSmth_twu*.nii']);

%% Define the data structure
dat=KORD_datapaths(rawdatapath,0);

%% Get all the ROIs from the data structure
ROIs_cell=[dat.ROIs.left dat.ROIs.right];
ROIstat=gk_getROIdata(statsFile,ROIs_cell);
ROI_names=fields(ROIstat)
ROIs_to_plot={ROI_names{5}, ROI_names{20}, ROI_names{12}, ROI_names{19}};

%% Iterate over subject files and get the data for each ROI
for i=1:numel(zScoreTcFilesHigh)
    zTCH_data(i)=gk_getROIdata(fullfile(zScoreTcFilesHigh(i).folder,zScoreTcFilesHigh(i).name),ROIs_cell);
end

for i=1:numel(zScoreTcFilesLow)
    zTCL_data(i)=gk_getROIdata(fullfile(zScoreTcFilesLow(i).folder,zScoreTcFilesLow(i).name),ROIs_cell);
end

%% Iterate over ROIs to plot and get the relevant data
% INFO: We check if the ROI is positive or negative on average and get the
% peak voxel in that direction

% change to true if you would like CSV file output
writeCSV=false;

for ri=1:numel(ROIs_to_plot)
    roiSign=sign(nanmean(ROIstat.(ROIs_to_plot{ri})));
    [~,maxind]=nanmax(roiSign*ROIstat.(ROIs_to_plot{ri}));
    for i=1:numel(zScoreTcFilesHigh)
        mTCH(i,:,ri)=zTCH_data(i).(ROIs_to_plot{ri})(maxind,:);
    end
    for i=1:numel(zScoreTcFilesLow)
        mTCL(i,:,ri)=zTCL_data(i).(ROIs_to_plot{ri})(maxind,:);
    end
    if writeCSV
        dlmwrite(fullfile(processedDataPath,[ROIs_to_plot{ri},'_High.csv']),mTCH(:,:,ri),'delimiter','\t')
        dlmwrite(fullfile(processedDataPath,[ROIs_to_plot{ri},'_Low.csv']),mTCL(:,:,ri),'delimiter','\t')
    end
end

%% Iterate over ROIs to plot and plot he mean plus SEM
figure; 
for ri=1:numel(ROIs_to_plot)
    subplot(2,2,ri); hold on;
    
    % LOW CONCENTRATION
    mL=mean(mTCL([1 3:12],:,ri));
    semL=std(mTCL([1 3:12],:,ri),1,1)/sqrt(size(mTCL([1 3:12],:,ri),1));
    errorbar(0.5:0.5:60,mL,semL,'m','LineWidth',0.5);
    plot(0.5:0.5:60,mean(mTCL([1 3:12],:,ri)),'m','LineWidth',3);
    % HIGH CONCENTRATION
    mH=mean(mTCH(:,:,ri));
    semH=std(mTCH(:,:,ri),1,1)/sqrt(size(mTCH,1));
    errorbar(0.5:0.5:60,mH,semH,'b','LineWidth',0.5);
    plot(0.5:0.5:60,mean(mTCH(:,:,ri)),'b','LineWidth',3);
    
    ylim([-20 20]);
end