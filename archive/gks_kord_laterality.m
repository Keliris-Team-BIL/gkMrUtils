% script to calculate the laterality indices


% Define the data paths
rawdatapath='/Users/gkeliris/Documents/DATA/Lore_KORD';
processedDataPath=cat(2,rawdatapath,'/ProcessedData_zscore');

% Define the data structure
dat=KORD_datapaths(rawdatapath,0);


% Define from which ROIs to get the data
ROIs_cell=[dat.ROIs.left dat.ROIs.right];

% Perform this once to create new AUC files with sign of stats
%kordFiles=dir([processedDataPath,'/KORD/HIGH_CONC/auc_z_adj_zscr_mSmth_twu*.nii']);
kordFiles=dir([processedDataPath,'/KORD/HIGH_CONC/auc_z_adj_zscr_mSmth_twu*.nii']);


%sgn_stats=sign(niftiread('~/Desktop/Zscore_smth3_95_80-100/Highconc_spmT_0001.nii'));
% % for i=1:numel(kordFiles)    
% %     data_info   = niftiinfo(fullfile(kordFiles(i).folder,kordFiles(i).name));
% %     data        = niftiread(data_info);
% % %     data1       = data.*sgn_stats;
% % %     niftiwrite(data1,fullfile(kordFiles(i).folder,['sgn_',kordFiles(i).name]),data_info);
% % %     data2       = -1*data;
% % %     niftiwrite(data2,fullfile(kordFiles(i).folder,['neg_',kordFiles(i).name]),data_info);

% % end

% Define the paths of the datafiles
%kordFiles=dir([processedDataPath,'/KORD/HIGH_CONC/neg_auc_z_adj_zscr_mSmth_twu*.nii']);
kordFiles=dir([processedDataPath,'/KORD/HIGH_CONC/auc_z_adj_zscr_mSmth_twu*.nii']);


% iterate each data file and get data
for i=1:numel(kordFiles)    
     ROI_data(i)=gk_getROIdata(fullfile(kordFiles(i).folder,kordFiles(i).name),ROIs_cell);
end

ROI_fields=fields(ROI_data);

left_idx=find(cell2mat(cellfun(@(x) ~isempty(strfind(x,'left')),ROI_fields,'UniformOutput',false)));
for i=1:numel(left_idx)
    ROI_name = ROI_fields{i}(1:end-5);
    L_ROI_name = [ROI_name '_left'];
    R_ROI_name = [ROI_name '_right'];
    for j=1:numel(ROI_data)
        Ldata=[ROI_data(:).(L_ROI_name)];
        Rdata=[ROI_data(:).(R_ROI_name)];
        res.(ROI_name).Ldata=Ldata;
        res.(ROI_name).Rdata=Rdata;
        
        for s=1:size(Ldata,2)
            clear Lpos Rpos sub_LI
            Lpos=Ldata(Ldata(:,s)>0, s);
            Rpos=Rdata(Rdata(:,s)>0, s);
            for li=1:numel(Lpos)
                Lt=sum(Lpos(Lpos>=Lpos(li)));
                Rt=sum(Rpos(Rpos>=Lpos(li)));
                sub_LI(li)=(Rt-Lt)/(Lt+Rt);
            end
            for ri=1:numel(Rpos)
                Lt=sum(Lpos(Lpos>=Rpos(ri)));
                Rt=sum(Rpos(Rpos>=Rpos(ri)));
                sub_LI(ri+numel(Lpos))=(Rt-Lt)/(Lt+Rt);
            end
            if exist('sub_LI','var')
                res.(ROI_name).aveLI(s)=mean(sub_LI);
            else
                res.(ROI_name).aveLI(s)=NaN;
            end
        end
       [res.(ROI_name).h res.(ROI_name).p]=ttest(res.(ROI_name).aveLI);
    end
     fprintf('AveLI\t=%.2f +- %.2f, p= %f (%s)\n',nanmean(res.(ROI_name).aveLI),nanstd(res.(ROI_name).aveLI)/sqrt(s),res.(ROI_name).p,ROI_name);
end
        
        
