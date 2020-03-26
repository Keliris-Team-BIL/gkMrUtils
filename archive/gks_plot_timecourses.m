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
        csvwrite(fullfile('~/Desktop/',[ROIs_to_plot{ri},'_High.csv']),mTCH(:,:,ri))
        csvwrite(fullfile('~/Desktop/',[ROIs_to_plot{ri},'_Low.csv']),mTCL(:,:,ri))
    end
end

%% Iterate over ROIs to plot and plot he mean plus SEM
figure; 
for ri=1:numel(ROIs_to_plot)
    subplot(2,2,ri); hold on;
    
    % LOW CONCENTRATION
    mL=mean(mTCL([1 3:12],:,ri));
    semL=std(mTCL([1 3:12],:,ri),1,1)/sqrt(size(mTCL([1 3:12],:,ri),1));
    errorbar(1:120,mL,semL,'m','LineWidth',0.5);
    plot(mean(mTCL([1 3:12],:,ri)),'m','LineWidth',3);
    % HIGH CONCENTRATION
    mH=mean(mTCH(:,:,ri));
    semH=std(mTCH(:,:,ri),1,1)/sqrt(size(mTCH,1));
    errorbar(1:120,mH,semH,'b','LineWidth',0.5);
    plot(mean(mTCH(:,:,ri)),'b','LineWidth',3);
    
    ylim([-20 20]);
end

