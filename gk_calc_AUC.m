function aucMap = gk_calc_AUC(dat,SHAMfolder,KORDfolder,R2_threshold,AUC_range,normtype)

%dat=KORD_datapaths(fileparts(pwd),0); % this works if you are in the code directory
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


% % PLOT THE AUC over a threshold 
% template=niftiread(dat.ROIs.mask{1});
% AUC_threshold=4;
% 
% % Plot each scan
% for sc=1:N_kord
%     aucSignif{sc}=goodFitIndex(find(abs(auc(:,sc)>AUC_threshold)));
%     mask=zeros(info.ImageSize);
%     mask(aucSignif{sc})=1;
%     figure;
%     for i=3:15 % the slices that are not empty
%         subplot(4,4,i);
%         image(100*cat(3,template(:,:,i)',template(:,:,i)',template(:,:,i)')); hold on;
%         im=imagesc(aucMap(:,:,i,sc)');
%         axis xy; caxis([-5 5])
%         im.AlphaData=mask(:,:,i)';
%     end
% end
% 
% % Plot the mean AUC of selected scans
% selectedScans=[1 3 4 6:12];
% av=mean(aucMap(:,:,:,selectedScans),4); 
% avSignif=zeros(info.ImageSize(1:3));
% avSignif(abs(av)>AUC_threshold)=1;
% figure;
% for i=3:15 % the slices that are not empty
%     subplot(4,4,i);
%     image(100*cat(3,template(:,:,i)',template(:,:,i)',template(:,:,i)')); hold on;
%     im=imagesc(av(:,:,i)');
%     axis xy; caxis([-5 5])
%     im.AlphaData=avSignif(:,:,i)';
% end
% print(fullfile(KORDfolder,'AUCmap'),'-dpdf','-fillpage')

