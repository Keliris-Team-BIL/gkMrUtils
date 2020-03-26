function gk_niftiM_to_percent(fileIn,fileOut,baseline)
% Usage: gk_niftiM_to_percent(fileIn,fileOut,baseline)
%
% fileIn    - the file to be loaded including full path
% fileOut   - the file to be saved including full path
% baseline  - the volumes to be used for the calculation of percent signal
%             (e.g [1:20])
%
% author: GAK

% load the data
info    = niftiinfo(fileIn);
data    = niftiread(fileIn);

%data=load_nii(fileIn);

% get the dimensions
dim=info.ImageSize;

% convert the data to 2D for easier calculations
data2D=reshape(data,prod(dim(1:3)),dim(4));

% calculate the mean of the baseline
s=mean(data2D(:,baseline)');
rep_s=repmat(s,dim(4),1);

% calculate percent based on baseline
percent=100*(((single(data2D')./rep_s))-1);
newimg=reshape(percent', dim);

% replace data in structure
%data.img=newimg;
%data.hdr.dime.datatype=16;
%data.hdr.dime.bitpix=32;

data = single(newimg);
info.Datatype='single';
info.BitsPerPixel=32;

% save to output file
%save_nii(data,fileOut);
niftiwrite(data,fileOut,info);
% give feedback
fprintf('Data saved to: %s\n',fileOut);