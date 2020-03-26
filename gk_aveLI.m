function res = gk_aveLI(niftifile_cell, left_ROIs_cell, right_ROIs_cell, verbose)
% USAGE: res = gk_aveLI(niftifile_cell, left_ROIs_cell, right_ROIs_cell, [verbose=true])
%
% INFO: This function takes as input a cell of nifti datafile paths and two cells
%       of Left and Right ROI nifti file paths to calculate the laterality using
%       the aveLI algorithm (Matsuo, K., Chen, S., Tseng, W. (2012). AveLI:
%       A robust lateralization index in functional magnetic resonance imaging
%       using unbiased threshold-free computation Journal of Neuroscience Methods
%       205(1), 119-129. https://dx.doi.org/10.1016/j.jneumeth.2011.12.020)
%       The index is calulated as (R-L)/(R+L) thus positive values indicate stronger
%       Right hemisphere activity
%
% INPUT:
% - niftifile_cell  : a cell of paths to nifti data files
% - left_ROIs_cell  : a cell of paths to nifti Left ROI files
% - right_ROIs_cell : a cell of paths to nifti Right ROI files
%   NOTE: the left and right ROI cells should have the same numbers of ROIs and
%         these ROIs should be ordered in the same way as they will be paired.
%         Further the ROI names are taken from the names of the ROI files. They
%         are assumed to be of the form ROIname_side_xxx or Side_ROIname_xx. The
%         side could be Left, Right, L, R and it is not case sensitive. The separator
%         can be underscore or space.
% - verbose         : if true (default) it will print the mean and sem of the AveLI
%                     over the datafiles for each ROI
%
% OUTPUT:
% - res   : a data structure that contains the aveLI per ROI and datafile
%
% AUTHOR: GAK, 12 Mar 2020


if nargin==3
  verbose=1;
end
if numel(left_ROIs_cell)~=numel(right_ROIs_cell)
  fprintf('ERROR: Number of left ROIs should be equal to number of right ROIs\n');
end

% iterate each data file and get data
for i=1:numel(niftifile_cell)
     L_ROI_data(i)=gk_getROIdata(niftifile_cell{i},left_ROIs_cell);
     R_ROI_data(i)=gk_getROIdata(niftifile_cell{i},right_ROIs_cell);
end
L_ROI_fields=fields(L_ROI_data);
R_ROI_fields=fields(R_ROI_data);


for i=1:numel(left_ROIs_cell)
    L_ROI_name = L_ROI_fields{i};
    R_ROI_name = R_ROI_fields{i};
    Lstr=regexp(upper(L_ROI_name),'(?<roi>\<[A-Z]*)[\s_](?<hem>LEFT[\>\s_]*|L[\>\s_]*|RIGHT[\>\s_]*|R[\>\s_]*).*|(?<hem>\<LEFT|\<L|\<RIGHT|\<R)[\s_](?<roi>[A-Z]*).*','names');
    Rstr=regexp(upper(R_ROI_name),'(?<roi>\<[A-Z]*)[\s_](?<hem>LEFT[\>\s_]*|L[\>\s_]*|RIGHT[\>\s_]*|R[\>\s_]*).*|(?<hem>\<LEFT|\<L|\<RIGHT|\<R)[\s_](?<roi>[A-Z]*).*','names');
    if ~strcmp(Lstr.roi,Rstr.roi)
      fprintf('Warning!!! Left and Right ROIs have different names (using name of LEFT)\n');
    end
    ROI_name=Lstr.roi;
    for j=1:numel(L_ROI_data)
        Ldata=[L_ROI_data(:).(L_ROI_name)];
        Rdata=[R_ROI_data(:).(R_ROI_name)];
        res.(Lstr.roi).Ldata=Ldata;
        res.(Lstr.roi).Rdata=Rdata;

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
    if verbose
     fprintf('AveLI = %+.2f +- %.2f, p = %f (%s)\n',nanmean(res.(ROI_name).aveLI),nanstd(res.(ROI_name).aveLI)/sqrt(s),res.(ROI_name).p,ROI_name);
    end
end
