function ROI_data = gk_getROIdata(data_nii,ROIs_cell)
% Usage: ROI_data = gk_getROIdata(nifti_filepath, ROIs_cell)
%
% INPUT:
% - nifti_filepath  : the full file path for the nifti file to get data from
% - ROIs_cell       : a cell of paths to ROI nifti files. For example you can
%                     join the left and right cells returned by the description
%                     file (i.e. ROIs_cell=[dat.ROIs.left dat.ROIs.right];)
%
% OUTPUT:
% - ROI_data        : a structure that contains the voxel data for each ROI
%
% author: GAK 3 Mar 2020


% load the data
data_info   = niftiinfo(data_nii);
data        = niftiread(data_nii);
dimData     = data_info.ImageSize;
if numel(dimData)==3; dimData(4)=1; end
data2D=reshape(data,prod(dimData(1:3)),dimData(4));
% change zero lines with NaNs
zero_lines=ismember(data2D,zeros(1,dimData(4)),'rows');
data2D(find(zero_lines),:)=NaN;


if ~iscell(ROIs_cell)
    if ischar(ROIs_cell)
        tmp=ROIs_cell;
        clear ROIs_cell;
        ROIs_cell{1}=tmp;
    else
        fprintf('ROIs_cell has to be a cell or a char containing paths to files\n')
        return
    end
end

for r=1:numel(ROIs_cell)
    [~,roi_name,~]=fileparts(ROIs_cell{r});
    ROI_info    = niftiinfo(ROIs_cell{r});
    ROI         = niftiread(ROIs_cell{r});

    % get the dimensions
    dimROI  = ROI_info.ImageSize;

    if isequal(dimData(1:3),dimROI)
        ROI_idx = find(ROI);
        ROI_data.(roi_name) = data2D(ROI_idx,:);
    end
end
