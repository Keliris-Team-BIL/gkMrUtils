function aucMap = gk_calc_AUC(dat,SHAMfolder,KORDfolder,R2_threshold,AUC_range,normtype)
% USAGE: aucMap = gk_calc_AUC(dat,SHAMfolder,KORDfolder,R2_threshold,AUC_range,normtype)
%
% INFO: This function calculates tha Area Under the Curve (AUC) values for a
%       specified range in the time course for the KORD experiments.
%
% INPUT:
% - dat         : the structure returned by KORD_datapaths
% - SHAMfolder  : the path to the folder with the SHAM data
% - KORDfolder  : the path to the folder with the KORD DATA
% - R2_threshold: franction (e.g. 0.95) only use the SHAM fits above this threshold
% - AUC_range   : the volumes to use for calculation of the AUC
% - normtype    : the pre created files to use (1 percent change, 2 zscore, 0 no norm)
%
% OUTPUT:
% - aucMap      : a variable with the AUC map. It also saves a .mat file.
% - NOTE: it will also create auc_XXX.nii files in the KORD directory per file.
%
% Author: GAK June 2019



if normtype==1
    adjstr='p_adj_';
elseif normtype==2
    adjstr='z_adj_';
else
    fprintf('normtype should be 1 if normalization was to percent signal or 2 if it was to z-score\n')
    return
end

% Load fitted SHAM and R2 map
info=niftiinfo(fullfile(SHAMfolder,'poly1_exp2_fit_mean_sham.nii'));
mean_sham_fit=niftiread(fullfile(SHAMfolder,'poly1_exp2_fit_mean_sham.nii'));
mean_sham_fit_R2=niftiread(fullfile(SHAMfolder,'R2_poly1_exp2_fit_mean_sham.nii'));
% select the voxels with fits> R2_threshold
goodFitIndex=find(mean_sham_fit_R2>R2_threshold);
% convert fitted sham to 2D for easier calculations and select good voxels
mean_sham_fit_2D=reshape(mean_sham_fit,prod(info.ImageSize(1:3)),info.ImageSize(4));
mean_sham_fit_2D_goodfit=mean_sham_fit_2D(goodFitIndex,:);

% Load all pre-processed kord scans
if normtype==1
    kordFiles=dir([KORDfolder 'pch_mSmth_twu*.nii']);
elseif normtype==2
    kordFiles=dir([KORDfolder 'zscr_mSmth_twu*.nii']);
end
N_kord=numel(kordFiles);
for i=1:N_kord
    kLow(:,:,:,:,i)=niftiread(fullfile(kordFiles(i).folder,kordFiles(i).name));
end
kLow_2Dall=reshape(kLow,prod(info.ImageSize(1:3)),info.ImageSize(4),N_kord);
kLow_2Dall_goodfit=kLow_2Dall(goodFitIndex,:,:);

% create the variable to store AUC values
aucMap=zeros([info.ImageSize(1:3) N_kord]);

% Iterate for each KORD scan
for sc=1:N_kord
    info=niftiinfo(fullfile(kordFiles(sc).folder,kordFiles(sc).name));
    if normtype==1 % percent
        % Regress each kord voxels baseline to the mean_sham_fit
        for k=1:numel(goodFitIndex)
            B(k,sc)=regress(kLow_2Dall_goodfit(k,1:20,sc)',mean_sham_fit_2D_goodfit(k,1:20)');
        end
        % Use the mean regression coefficient to subtract sham for each KORD
        mB=mean(B,1);
    elseif normtype==2
        mB=ones(N_kord,1);
    end
    mB_kLow_2Dall_goodfit(:,:,sc)=kLow_2Dall_goodfit(:,:,sc)-(mean_sham_fit_2D_goodfit.*repmat(mB(sc),1,120));

    % Save niftis
    kord=zeros(info.ImageSize);
    kLow_2Dall(goodFitIndex,:,sc)=mB_kLow_2Dall_goodfit(:,:,sc);
    kord=reshape(kLow_2Dall(:,:,sc),info.ImageSize);
    niftiwrite(kord,fullfile(kordFiles(sc).folder,[adjstr, kordFiles(sc).name]),info);

    % Calculate area under the curve (AUC)
    temp=zeros(info.ImageSize(1:3));
    auc(:,sc)=mean(mB_kLow_2Dall_goodfit(:,AUC_range,sc),2);
    temp(goodFitIndex)=auc(:,sc);
    aucMap(:,:,:,sc)=temp;
    info=niftiinfo(dat.ROIs.mask{1});
    info.Datatype='single';
    info.BitsPerPixel=16;

    % Save niftis
    niftiwrite(single(aucMap(:,:,:,sc)),fullfile(kordFiles(sc).folder,['auc_',adjstr, kordFiles(sc).name]),info);


end
save(fullfile(KORDfolder,'AUCmap'),'aucMap');
