function gk_kord_preprocess(filenames,outfolder,mask,smth,bsl,normtype)
% Usage: gk_kord_preprocess(filenames,outfolder,mask,smth,bsl,normtype)
%
% function to preprocess the data of the KORD project
%
% INPUT
% filenames:    a cell structure with the paths to the files e.g. the one
%               stored in the dat returned by the description file
%               (KORD_datapaths) e.g: dat.sham.lowconc.normalized
%
% outfolder:    the folder where the preprocessed files will be written
%
% mask:         the location of the mask nifti e.g the one stored in the
%               description file (KORD_datapaths) e.g. dat.ROIs.mask{1}
%
% smth:         the smoothing kernel, can be scalar for isometric or a
%               vector of 3 [sx, sy, sz] for 3D. leave empty [] for none
%
% bsl:          the baseline volumes for calculating percent signal change,
%               leave empty [] for none
%
% normtype:     1 convert to percent change, 2 convert to zscore, 0 no norm
%
% example:
%   dat=KORD_datapaths;
%   processedDataPath=cat(2,fileparts(pwd),'/Test');
%   vox_x=0.3;
%   gk_kord_preprocess(dat.sham.lowconc.normalized,processedDataPath,dat.ROIs.mask{1},2*vox_x,1:20);
%
% GAK May 2019




for s=1:numel(filenames)    
    fileIn=filenames{s};
    [~,fname,fextension] = fileparts(fileIn);
    
    % SMOOTHING
    if ~isempty(smth)
        smOut=fullfile(outfolder,cat(2,'Smth_',fname,fextension));
        spm_smooth(fileIn,smOut,smth);
        fname=cat(2,'Smth_',fname);
    end
    % MASKING
    if ~isempty(mask)
        spm_mask(mask,fullfile(outfolder,cat(2,fname,fextension)));
        fname=cat(2,'m',fname);
    end
    % Normalize to PERCENT CHANGE or ZSCORE RELATIVE TO BASELINE (bsl)
    if ~isempty(bsl) && normtype==1
        fileOut=fullfile(outfolder,cat(2,'pch_',fname,fextension));
        gk_niftiM_to_percent(fullfile(outfolder,cat(2,fname,fextension)),fileOut,bsl);
    elseif ~isempty(bsl) && normtype==2
        fileOut=fullfile(outfolder,cat(2,'zscr_',fname,fextension));
        gk_niftiM_to_zscore(fullfile(outfolder,cat(2,fname,fextension)),fileOut,bsl);
    end
end
