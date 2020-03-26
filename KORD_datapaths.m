function dat = KORD_datapaths(rawDataPath,verbose)
% USAGE: dat = KORD_datapaths(rawDataPath)
%
% rawDataPath : full path to the raw data folder without / at the end

% Note: always use forward slashes (/) instead of backslashes (\) in path
% definitions as they work in both windows and linux platforms.

if ~exist('verbose','var')
    verbose=1;
end

%%% PATHS THAT POTENTIALLY NEED TO BE UPDATED %%%%%%%%%%%%%%%%%%%%%%%%%%%%
pathLowConc         =cat(2,rawDataPath,'/PhMRI_LowConc/');
pathHighConc        =cat(2,rawDataPath,'/PhMRI_HighConc/');
pathLowConcPerc     =cat(2,pathLowConc,'Percentchange/');
pathHighConcPerc    =cat(2,pathHighConc,'Pecentchange/');
pathROIs            =cat(2,rawDataPath,'/ROIs/');


% ANIMAL IDs
dat.sham.IDs={'1284_LL','1284_R','1284_RL','1284_RR','1286_L','1286_R'};
dat.kord.IDs={'1288_L','1288_LL','1288_NM','1288_R','1288_RL','1288_RR',...
             '1290_L','1290_LL','1290_NM','1290_R','1290_RL','1290_RR'};
         
% 
%% RAW DATA
% Raw data SHAMs: Low concentration
dat.sham.lowconc.raw = {...
    cat(2,pathLowConc,dat.sham.IDs{1},'/scan9/1284_LL9.nii'),...
    cat(2,pathLowConc,dat.sham.IDs{2},'/scan9/1284_R9.nii'),...
    cat(2,pathLowConc,dat.sham.IDs{3},'/scan10/1284_RL10.nii'),...
    cat(2,pathLowConc,dat.sham.IDs{4},'/scan8/1284_RR8.nii'),...
    cat(2,pathLowConc,dat.sham.IDs{5},'/scan8/1286_L8.nii'),...
    cat(2,pathLowConc,dat.sham.IDs{6},'/scan10/1286_R10.nii'),...
    };

% Raw data KORDs: Low concentration
dat.kord.lowconc.raw = {...
    cat(2,pathLowConc,dat.kord.IDs{1},'/scan9/1288_L9.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{2},'/scan8/1288_LL8.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{3},'/scan10/1288_NM10.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{4},'/scan8/1288_R8.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{5},'/scan8/1288_RL8.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{6},'/scan9/1288_RR9.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{7},'/scan11/1290_L11.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{8},'/scan10/1290_LL10.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{9},'/scan8/1290_NM8.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{10},'/scan12/1290_R12.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{11},'/scan13/1290_RL13.nii')...
    cat(2,pathLowConc,dat.kord.IDs{12},'/scan8/1290_RR8.nii')...
    };

% Raw data KORDs: High concentration
dat.kord.highconc.raw = {...
    cat(2,pathHighConc,dat.kord.IDs{1},'/scan8/1288_L8.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{2},'/scan9/1288_LL9.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{3},'/scan8/1288_nm8.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{4},'/scan8/1288_R8.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{5},'/scan8/1288_RL8.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{6},'/scan9/1288_RR9.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{7},'/scan10/1290_L10.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{8},'/scan8/1290_LL8.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{9},'/scan8/1290_nm8.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{10},'/scan8/1290_R8.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{11},'/scan9/1290_RL9.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{12},'/scan8/1290_RR8.nii')...
    };

%% NORMALIZED DATA
% Normalized data SHAMs: Low concentration
dat.sham.lowconc.normalized = {...
    cat(2,pathLowConc,dat.sham.IDs{1},'/Normalized_scan/twu1284_LL9.nii'),...
    cat(2,pathLowConc,dat.sham.IDs{2},'/Normalized_scan/twu1284_R9.nii'),...
    cat(2,pathLowConc,dat.sham.IDs{3},'/Normalized_scan/twu1284_RL10.nii'),...
    cat(2,pathLowConc,dat.sham.IDs{4},'/Normalized_scan/twu1284_RR8.nii'),...
    cat(2,pathLowConc,dat.sham.IDs{5},'/Normalized_scan/twu1286_L8.nii'),...
    cat(2,pathLowConc,dat.sham.IDs{6},'/Normalized_scan/twu1286_R10.nii'),...
    };

% Normalized data KORDs: Low concentration
dat.kord.lowconc.normalized = {...
    cat(2,pathLowConc,dat.kord.IDs{1},'/Normalized_scan/twu1288_L9.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{2},'/Normalized_scan/twu1288_LL8.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{3},'/Normalized_scan/twu1288_nm10.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{4},'/Normalized_scan/twu1288_R8.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{5},'/Normalized_scan/twu1288_RL8.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{6},'/Normalized_scan/twu1288_RR9.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{7},'/Normalized_scan/twu1290_L11.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{8},'/Normalized_scan/twu1290_LL10.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{9},'/Normalized_scan/twu1290_NM8.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{10},'/Normalized_scan/twu1290_R12.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{11},'/Normalized_scan/twu1290_RL13.nii')...
    cat(2,pathLowConc,dat.kord.IDs{12},'/Normalized_scan/twu1290_RR8.nii')...
    };

% Normalized data KORDs: High concentration
dat.kord.highconc.normalized={...
    cat(2,pathHighConc,dat.kord.IDs{1},'/Normalized_scan/twu1288_L8.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{2},'/Normalized_scan/twu1288_LL9.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{3},'/Normalized_scan/twu1288_nm8.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{4},'/Normalized_scan/twu1288_R8.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{5},'/Normalized_scan/twu1288_RL8.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{6},'/Normalized_scan/twu1288_RR9.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{7},'/Normalized_scan/twu1290_L10.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{8},'/Normalized_scan/twu1290_LL8.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{9},'/Normalized_scan/twu1290_nm8.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{10},'/Normalized_scan/twu1290_R8.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{11},'/Normalized_scan/twu1290_RL9.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{12},'/Normalized_scan/twu1290_RR8.nii')...
    };

%% NORMALIZED DATA - MASKED
% Normalized data SHAMs: Low concentration
dat.sham.lowconc.norm_masked = {...
    cat(2,pathLowConc,dat.sham.IDs{1},'/Normalized_scan/twu1284_LL_masked.nii'),...
    cat(2,pathLowConc,dat.sham.IDs{2},'/Normalized_scan/twu1284_R_masked.nii'),...
    cat(2,pathLowConc,dat.sham.IDs{3},'/Normalized_scan/twu1284_RL_masked.nii'),...
    cat(2,pathLowConc,dat.sham.IDs{4},'/Normalized_scan/twu1284_RR_masked.nii'),...
    cat(2,pathLowConc,dat.sham.IDs{5},'/Normalized_scan/twu1286_L_masked.nii'),...
    cat(2,pathLowConc,dat.sham.IDs{6},'/Normalized_scan/twu1286_R_masked.nii'),...
    };

% Normalized data KORDs: Low concentration
dat.kord.lowconc.norm_masked = {...
    cat(2,pathLowConc,dat.kord.IDs{1},'/Normalized_scan/twu1288_L_masked.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{2},'/Normalized_scan/twu1288_LL_masked.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{3},'/Normalized_scan/twu1288_NM_masked.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{4},'/Normalized_scan/twu1288_R_masked.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{5},'/Normalized_scan/twu1288_RL_masked.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{6},'/Normalized_scan/twu1288_RR_masked.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{7},'/Normalized_scan/twu1290_L_masked.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{8},'/Normalized_scan/twu1290_LL_masked.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{9},'/Normalized_scan/twu1290_NM_masked.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{10},'/Normalized_scan/twu1290_R_masked.nii'),...
    cat(2,pathLowConc,dat.kord.IDs{11},'/Normalized_scan/twu1290_RL_masked.nii')...
    cat(2,pathLowConc,dat.kord.IDs{12},'/Normalized_scan/twu1290_RR_masked.nii')...
    };

% Normalized data KORDs: High concentration
dat.kord.highconc.norm_masked={...
    cat(2,pathHighConc,dat.kord.IDs{1},'/Normalized_scan/twu1288_L_masked.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{2},'/Normalized_scan/twu1288_LL_masked.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{3},'/Normalized_scan/twu1288_NM_masked.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{4},'/Normalized_scan/twu1288_R_masked.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{5},'/Normalized_scan/twu1288_RL_masked.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{6},'/Normalized_scan/twu1288_RR_masked.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{7},'/Normalized_scan/twu1290_L_masked.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{8},'/Normalized_scan/twu1290_LL_masked.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{9},'/Normalized_scan/twu1290_NM_masked.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{10},'/Normalized_scan/twu1290_R_masked.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{11},'/Normalized_scan/twu1290_RL_masked.nii'),...
    cat(2,pathHighConc,dat.kord.IDs{12},'/Normalized_scan/twu1290_RR_masked.nii')...
    };

% Manually drawn ROIs
dat.ROIs.left={...
    cat(2,pathROIs,'ACA_left_GAK.nii'),...
    cat(2,pathROIs,'Auditoryctx_left.nii'),...
    cat(2,pathROIs,'Hippocampus_left.nii'),...
    cat(2,pathROIs,'Insular_left.nii'),...
    cat(2,pathROIs,'SScortex_left.nii'),...
    cat(2,pathROIs,'Striatum_left.nii'), ...
    cat(2,pathROIs,'Thalamus_left.nii'),...
    cat(2,pathROIs,'Visualcortex_left.nii'),...
    cat(2,pathROIs,'PIR_left.nii'),...
    cat(2,pathROIs,'RSP_left.nii'),...
    cat(2,pathROIs,'vPAL_left.nii'),...
    cat(2,pathROIs,'BLA_Left_GAK.nii')...
    };

%cat(2,pathROIs,'ACA_right.nii'),...
dat.ROIs.right={
    cat(2,pathROIs,'ACA_right.nii'),...
    cat(2,pathROIs,'Auditoryctx_right.nii'),...
    cat(2,pathROIs,'Hippocampus_right.nii'),...
    cat(2,pathROIs,'Insular_right.nii'),...
    cat(2,pathROIs,'SScortex_right.nii'),...
    cat(2,pathROIs,'Striatum_right.nii'),...
    cat(2,pathROIs,'Thalamus_right_GAK.nii'),...
    cat(2,pathROIs,'Visualcortex_right.nii'),...
    cat(2,pathROIs,'ACA_Right_january.nii'),...
    cat(2,pathROIs,'ACA_Right_january_3.nii'),...
    cat(2,pathROIs,'PIR_right.nii'),...
    cat(2,pathROIs,'RSP_right.nii'),...
    cat(2,pathROIs,'vPAL_right.nii')...
    };

dat.ROIs.mask={cat(2,pathROIs,'brainmask.nii')};


%% VERIFY THAT ALL FILES IN STRUCTURE ARE EXISTING
% Note: exclude fields that are not files (e.g. IDs)
if verbose
    clear global;
    global s;
    unfold(dat,false);
    
    for si=1:numel(s)
        %    strfind(s{si},'IDs')
        if isempty(strfind(s{si},'IDs'))
            for i=1:numel(eval(s{si}))
                if ~isfile(eval([s{si},'{',num2str(i),'}']))
                    fprintf('\nWARNING! NOT FOUND File: %s\n\n',...
                        eval([s{si},'{',num2str(i),'}']));
                else
                    fprintf('VERIFIED File: %s\n',eval([s{si},'{',num2str(i),'}']));
                end
            end
        end
    end
    clear global
end
