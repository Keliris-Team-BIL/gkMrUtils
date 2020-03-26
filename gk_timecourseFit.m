function ft = gk_timecourseFit(fileIn,maskIn,plotfit,params)
% Usage: ft = gk_timecourseFit(fileIn,maskIn,[plotfit],[params])
%
% params    -struct with fields: 
%               .baseline (the baseline period)              def: [1:20]
%               .x1 (the points used for the linear fit)     def: [5:20]
%               .x2 (the points used for the double exp fit) def: [25:120]
%               .x3 (where the crossing is expected)         def: [21:30]
% 
% GAK May 2019
if ~exist('plotfit','var') || isempty(plotfit)
    plotfit=false;
end

if ~exist('params','var')
    bsl=1:20;
    x1=5:20;
    x2=25:120;
    x3=21:30;
else
    bsl=params.baseline;
    x1=params.x1;
    x2=params.x2;
    x3=params.x3;
end


[filepath,fname,fextension] = fileparts(fileIn);
% load scan
sInfo=niftiinfo(fileIn);
nTimepoints=sInfo.ImageSize(4);
s=niftiread(fileIn);
dime=sInfo.ImageSize;
% load mask
mInfo=niftiinfo(maskIn);
mask=niftiread(maskIn);
maskind=find(mask);

% convert scan to 2D for easier calculations
s2D=reshape(s,prod(dime(1:3)),dime(4));
s2Dfit=single(zeros(size(s2D)));

% Use a linear fit for the baseline
ft1 = fittype('poly1');
% Use a double exponential for after the injection
ft2 = fittype('exp2');
% the options below are to force the first exponential to be negative
fo2 = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[-Inf,-Inf,-Inf,-Inf],...
               'Upper',[Inf,0,Inf,Inf],...
               'StartPoint',[0,-1,0,0]);


%% iterate over the voxels in mask and fit
Nvoxels=numel(maskind);
h=waitbar(0,'Processing...');
tic
for i=1:Nvoxels
    msg=['Processing voxel ',num2str(i),'/',num2str(Nvoxels)];
    waitbar(i/Nvoxels,h,msg)
    % estimate fitting params
    [ft(i).f1,ft(i).S1]=fit(x1',double(s2D(maskind(i),x1))',ft1);
    [ft(i).f2,ft(i).S2]=fit(x2',double(s2D(maskind(i),x2))',ft2,fo2);
    
    % find where the two fitted lines cross
    cross_ind=find(ft(i).f1(x3) > ft(i).f2(x3));
    if isempty(cross_ind)
        cross_ind=find(ft(i).f1(x3) < ft(i).f2(x3));
        cross_ind=max(cross_ind);
    end
    % apply the fit to estimate curve
    s2Dfit(maskind(i),:)=[ft(i).f1(1:bsl(end)+cross_ind(1)); ft(i).f2(x3(1)+cross_ind(1):nTimepoints)]';
end
toc
close(h);
% prepend and save nifti file
niftiwrite(reshape(s2Dfit,dime),fullfile(filepath,['poly1_exp2_fit_',fname,fextension]),sInfo);
save(fullfile(filepath,['poly1_exp2_fit_',fname]),'ft');
% Get adjusted R^2 stats
temp=[ft.S2];
exp2_R2=[temp.adjrsquare];

% create mask with goodness of fit values for easier selection later
mask=single(mask);
mask(maskind)=exp2_R2;
mInfo.Datatype='single';
mInfo.BitsPerPixel=32;
mInfo.Description='Rsquare_of_exp2_fit';
niftiwrite(mask,fullfile(filepath,['R2_poly1_exp2_fit_',fname,fextension]),mInfo);


% Plot PDF files of the fit if plotfit is true
if plotfit
    % sort the voxels from best to worst
    [~, sortedvox]=sort(exp2_R2,'descend');
    % Plot all based on R2
    figure('Position',[0 0 700 1000])
    for v=1:Nvoxels
        if rem(v,70); sp=rem(v,70); else; sp=70; end
        subplot(10,7,sp); hold on;
        plot(s2D(maskind(sortedvox(v)),:)','.');
        plot(s2Dfit(maskind(sortedvox(v)),:)');
        if sp==70
            print(num2str(ceil(v/70)),'-dpdf','-fillpage')
            clf;
        end
    end
    if sp<70
        print(num2str(ceil(v/70)),'-dpdf','-fillpage')
    end
    close(gcf)
end

