function aveli()

% Welcome to AveLI!
% This script computes a lateralization index called AveLI. 
% In this script, you will specify 3 types of files: data file, left and right ROI mask files. 
% Please be aware that you have to prepare these 3 files in the same brain template space
% such as the Montreal Neurological Institute (MNI) space.
% 
% In advance, you need to conduct the first level analysis by SPM software first.
% You will then obtain spmT_*.img as well as con_*.img files according to the contrast you make. 
% 
% You can analyze multiple sessions/subjects. If it is the case, all data files should have 
% the same dimensions in the same template space (e.g., 91*109*91 in the MNI space).
% 
% The data file you specify in this script may be the spmT_*.img, which contains t-values 
% in the voxels. 
% You may instead choose the con_*.img, which contains beta values (i.e., contrast estimates)
% in the voxels. 
% Other files would be also possible. Please select what you consider the best suited. 
% 
% For the ROI mask files, you can choose binary files (i.e., contain either 0 or 1 in the 
% voxels). The voxel value 1 indicates the locations of anatomy structures you are interested  
% in (e.g., inferior frontal gyrus). The lateralization index will be computed only using the 
% voxels with the value 1. 
% 
% The files you prepare as the ROI files may not be binary files and may have various values 
% in the voxels other than 0 or 1.
% Also, the files may have different resolutions from ones in the data file.  
% Even in these cases, this script will make new binary files that have the same resolutions  
% with the data files by applying a cut-off threshold at 0.5. 
% Please modify the threshold if necessary. 
% 
% We appreciate for you to cite the following paper: 
% Matsuo K, Chen SH, Tseng WY. 2012. AveLI: A robust lateralization index in functional magnetic resonance imaging using unbiased threshold-free computation. J Neurosci Methods. 205(1):119-129.
% The pubmed website is as follows: http://www.ncbi.nlm.nih.gov/pubmed/22233778
% 
% The download website of the AveLI script is as follows: 
% http://aveli.web.fc2.com/
% 
% If you find any bugs and problems, please let me know.
% Kayako Matsuo, Ph.D.
% Email:  info.aveli [at] gmail.com
% 
% Thank you very much for your interest in AveLI.
% 
% 
% History
% 2012.3.21. Uploaded first by Kayako @ NTUH, Taiwan, written on March 20, 2012
% 2015.1.29. Uploaded modified version by Kayako @ CPMMH, HUSM, Hamamatsu, Japan. 
% 2016.10.27. Modified by Kayako @ DMU, Tochigi, Japan
% 2017.4.3.  Modified (mask spec) by Kayako @ DMU, Tochigi, Japan
% 2019.12.6.  Comment out spm_imcalc_ui by Kayako @ DMU, Tochigi, Japan

tic

% ******************* Preparations *********************
% check for SPM 
  if isempty(which('spm')),
    error('No "spm" found, please addpath an SPM folder');
  end

% ****** File calling ******
% Data files
  fff=spm_select(Inf,'image','Select data file(s) to compute LI (multiple sessions OK)');

  % no file specification causes error
    if isempty(fff)
      error('Please select data file(s)!')
    end

  % check file number of fff
    loopn_fff=size(fff,1);

  % obtain dimensions of fff
    fffdim = zeros(loopn_fff,3);
    for i = 1:loopn_fff
      fnow = spm_vol(fff(i,:));
      fffdim(i,:) = fnow.dim;
    end

  % check equivalence
    if loopn_fff==1
      fffdimave = fffdim;

    else
    fffdimave = mean(fffdim);
    for i = 1:loopn_fff
      if ~isequal(fffdimave, fffdim(i,:))
        error('Sorry, the voxel dimensions should be the same for all data files!')
      end
    end

    end

% Left ROI file
  leftmaskf=spm_select(1,'image','Select a file for the Left ROI');

  % no file specification causes error
    if isempty(leftmaskf)
      error('Please select a Left ROI file!');
    end

% Right ROI file
  rightmaskf=spm_select(1,'image','Select a file for the Right ROI');

  % ask to choose different files for left and right
    while isequal(rightmaskf, leftmaskf)
      rightmaskf=spm_select(1,'image','Right should be different from Left!');
    end

  % no file specification causes error
    if isempty(rightmaskf)
      error('Please select a Right ROI file!');
    end


% ****** Create new ROI files to use ******
% define file name for ROIs to be saved
  % Left mask 
    [pth_l, nam_l, ext_l] = fileparts(leftmaskf);
    
    % ext
      if isequal(ext_l,'.nii,1')
        ext_l = '.nii';
      elseif isequal(ext_l,'.img,1')
        ext_l = '.img';
      end

    leftmask = strcat(nam_l,'_aveli_used',ext_l);

  % Right mask
    [pth_r, nam_r, ext_r] = fileparts(rightmaskf);
    
    % ext
      if isequal(ext_r,'.nii,1')
        ext_r = '.nii';
      elseif isequal(ext_r,'.img,1')
        ext_r = '.img';
      end

    rightmask = strcat(nam_r,'_aveli_used',ext_r);

% temporary nulvol to make the resolutions same
  % make nulvol
    nulvol = spm_vol(fff(1,:));
    nuls = zeros(fffdimave);

  % Left
    leftnulfname = strcat('nulvol_left',ext_l);
    nulvol.fname = leftnulfname;
    spm_write_vol(nulvol, nuls);

  % Right
    rightnulfname = strcat('nulvol_right',ext_r);
    nulvol.fname = rightnulfname;
    spm_write_vol(nulvol, nuls);

% Make combination of files for ImCalc 
  combi_l = {leftnulfname;deblank(char(leftmaskf))};
  combi_r = {rightnulfname;deblank(char(rightmaskf))};

% Execute ImCalc
  % I will comment out spm_imcalc_ui (2019.12.6.)
  % I will leave older comment (the following 2 lines)
    % if you have spm_imcalc_ui (as in SPM8), it is used
    % if you only have spm_imcalc (as in SPM12), then it is used

%  askui=which('spm_imcalc_ui');

%  if isempty(askui)
  leftmask = spm_imcalc(combi_l, leftmask, 'i1-i1+(i2>0.5)');
  rightmask = spm_imcalc(combi_r, rightmask, 'i1-i1+(i2>0.5)');

%  else
%  leftmask = spm_imcalc_ui(combi_l, leftmask, 'i1-i1+(i2>0.5)',{0,0,typeoffff1,1});
%  righttmask = spm_imcalc_ui(combi_r, rightmask, 'i1-i1+(i2>0.5)',{0,0,typeoffff1,1});

%  end


% delete nulvol
  delete nulvol*

% Take indices of Left and Right
  ind_left = find(spm_read_vols(spm_vol(leftmask.fname)));
  ind_right = find(spm_read_vols(spm_vol(rightmask.fname)));


% ****** Preparation of the summary ******
% By each session/subject
% number of voxels with positive values
  Sum4VX_left=zeros(loopn_fff,1);
  Sum4VX_right=zeros(loopn_fff,1);

% summation of voxel intensity values
  Sum4value_left=zeros(loopn_fff,1);
  Sum4value_right=zeros(loopn_fff,1);

% average voxel values
  Sum4AveVoxel_left=zeros(loopn_fff,1);
  Sum4AveVoxel_right=zeros(loopn_fff,1);

% lateralization index by voxel values
  Sum4AveLI_value=zeros(loopn_fff,1);
  Sum4baseLI_value=zeros(loopn_fff,1);

% lateralization index by voxel number
  Sum4AveLI_vn=zeros(loopn_fff,1);
  Sum4baseLI_vn=zeros(loopn_fff,1);


% ****** Preparation of the text log file *****
% define file name 
  cc=clock; ddd=date; 
  yy=cc(1); yystr=num2str(yy);
  mm=cc(2); mmstr=num2str(mm);
  dd=cc(3); ddstr=num2str(dd);
  hh=cc(4); hhstr=num2str(hh);
   if length(hhstr)<2 hhstr=strcat('0',hhstr); end;
  mm=cc(5); mmstr=num2str(mm);
   if length(mmstr)<2 mmstr=strcat('0',mmstr); end;
  ss=cc(6); ssstr=num2str(ss);
  fnamelog=strcat('AveLI-',ddd,'-',hhstr,mmstr,'.',ssstr,'.txt');

  fid = fopen(fnamelog,'w');


% Header of the log 
  fprintf(fid,'***** Welcome to AveLI *****');fprintf(fid,'\n');

% Information of specified files
  fprintf(fid, 'Information of files you specified');fprintf(fid,'\n');
  fprintf(fid, 'ROI files');fprintf(fid,'\n');
  fprintf(fid, 'Left ROI');fprintf(fid,'\t');
  fprintf(fid, '%s', leftmaskf);fprintf(fid,'\n');
  fprintf(fid, 'Right ROI');fprintf(fid,'\t');
  fprintf(fid, '%s', rightmaskf);fprintf(fid,'\n');
  fprintf(fid, 'Data files');fprintf(fid,'\n');

  for i = 1:loopn_fff
    ssnum = num2str(i);
    fprintf(fid, '%s', ssnum);fprintf(fid,'\t');
    fprintf(fid, '%s', fff(i,:));fprintf(fid,'\n');
  end   % end of i

  fprintf(fid,'\n');

% Header of Results table
% order: SS, AveLI, AveLI(voxel), baseLI, baseLI(voxel), Lt voxels, Rt voxels, 
%        Lt sum values, Rt sum values, Lt mean value, Rt mean value
  fprintf(fid, 'Results');fprintf(fid,'\n');
  fprintf(fid, 'SS');fprintf(fid,'\t');
  fprintf(fid, 'AveLI');fprintf(fid,'\t');
  fprintf(fid, 'AveLI(voxel)');fprintf(fid,'\t');
  fprintf(fid, 'baseLI');fprintf(fid,'\t');
  fprintf(fid, 'baseLI(voxel)');fprintf(fid,'\t');
  fprintf(fid, 'Lt voxels');fprintf(fid,'\t');
  fprintf(fid, 'Rt voxels');fprintf(fid,'\t');
  fprintf(fid, 'Lt sum values');fprintf(fid,'\t');
  fprintf(fid, 'Rt sum values');fprintf(fid,'\t');
  fprintf(fid, 'Lt mean value');fprintf(fid,'\t');
  fprintf(fid, 'Rt mean value');fprintf(fid,'\n');

% display message
  disp(' ');
  disp('Files prepared.');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ************** Main Loop Start **************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  for sslp = 1:loopn_fff   % this is for ss

  % call data for now
    vnow = spm_read_vols(spm_vol(deblank(char(fff(sslp,:)))));

  % change NaN to zero
    nanloc = isnan(vnow);
    vnow(nanloc) = 0;

  % change negative to zero
    negloc = vnow < 0;
    vnow(negloc) = 0;

  % extract left and right ROI areas, only positive values
    left_now = vnow(ind_left);
    right_now = vnow(ind_right);

    left_posi = find(left_now>0);
    right_posi = find(right_now>0);

    left_now = left_now(left_posi);
    right_now = right_now(right_posi);

    left_now = flipud(sort(left_now));
    right_now = flipud(sort(right_now));

  % voxel number 
    Sum4VX_left(sslp) = size(left_now,1);
    Sum4VX_right(sslp) = size(right_now,1);

  % sum of voxel intensities
    Sum4value_left(sslp) = sum(left_now);
    Sum4value_right(sslp) = sum(right_now);

  % average voxel values
    Sum4AveVoxel_left(sslp) = mean(left_now);
    Sum4AveVoxel_right(sslp) = mean(right_now);


% ***************************************************************
% ****************** #1 Calculation of AveLI ********************
% ***************************************************************

  % Indicate the start of the computation
    disp('Now computing... ');

  % to find the temporary threshold
    both_now = [left_now;right_now];
    both_now = flipud(sort(both_now));

  % mind no activation in the ROIs
    if isempty(both_now)    
      left_value_tmp = 0;
      left_vx_tmp = 0;
      right_value_tmp = 0;
      right_vx_tmp = 0;
      wlit = NaN;
      wliv = NaN;

    else

  % initialize spm progress bar
    spm_progress_bar('Init',100, 'each LI computation');

  % count loop number for voxels within both ROIs
    allvxnum = size(both_now,1);

  % prepare wlit for voxel value LI, wliv for voxel number LI
    wlit = zeros(allvxnum,1);
    wliv = zeros(allvxnum,1);

  % Start of calculation of each LI
    for eachli = 1:allvxnum

      thres_now = both_now(eachli);

    % initialization
      left_value_tmp = 0;
      left_vx_tmp = 0;
      right_value_tmp = 0;
      right_vx_tmp = 0;


    % temporary values for Left and Right
    % Left
      left_value_tmp = left_now(left_now >= thres_now);
      if isempty(left_value_tmp)
        left_value_tmp = 0;
        left_vx_tmp = 0;

      else
        left_vx_tmp = size(left_value_tmp,1);
        left_value_tmp = sum(left_value_tmp);

      end

    % Right
      right_value_tmp = right_now(right_now >= thres_now);
      if isempty(right_value_tmp)
        right_value_tmp = 0;
        right_vx_tmp = 0;

      else
        right_vx_tmp = size(right_value_tmp,1);
        right_value_tmp = sum(right_value_tmp);

      end


    % temporary LI (each li)
      wlit(eachli) = (left_value_tmp - right_value_tmp) / (left_value_tmp + right_value_tmp);
      wliv(eachli) = (left_vx_tmp - right_vx_tmp) / (left_vx_tmp + right_vx_tmp);


    % for spm progress bar
      epv100=round((eachli/allvxnum)*100);

      if epv100 >= 100
          ppp=100;
        elseif epv100 >= 90
          ppp=90;
        elseif epv100 >= 80
          ppp=80;
        elseif epv100 >= 70
          ppp=70;
        elseif epv100 >= 60
          ppp=60;
        elseif epv100 >= 50
          ppp=50;
        elseif epv100 >= 40
          ppp=40;
        elseif epv100 >= 30
          ppp=30;
        elseif epv100 >= 20
          ppp=20;
        elseif epv100 >= 10
          ppp=10;
        else
          ppp=0;
      end;


      spm_progress_bar('Set',ppp);


    end    % end of eachli

      spm_progress_bar('Clear');


    end      % end of mind no activation in the ROIs

  % Computation of AveLI 
    Sum4AveLI_value(sslp) = mean(wlit);
    Sum4AveLI_vn(sslp) = mean(wliv);



% ********************************************
% ************** #2 Raw LI (baseLI) **********
% ********************************************

  % No weighting LI using all voxels (t/beta > 0), that is, baseLI, for your infomation

    Sum4baseLI_value(sslp) = (left_value_tmp - right_value_tmp) / (left_value_tmp + right_value_tmp);
    Sum4baseLI_vn(sslp) = (left_vx_tmp - right_vx_tmp) / (left_vx_tmp + right_vx_tmp);



% ******************************************************
% Writing summary
% ******************************************************

  % order: SS, AveLI, AveLI(voxel), baseLI, baseLI(voxel), Lt voxels, Rt voxels, 
  %        Lt sum values, Rt sum values, Lt mean value, Rt mean value

    ssnum = num2str(sslp);
    fprintf(fid, '%s', ssnum);fprintf(fid, '\t');
    fprintf(fid, '%f', Sum4AveLI_value(sslp));fprintf(fid, '\t');
    fprintf(fid, '%f', Sum4AveLI_vn(sslp));fprintf(fid, '\t');
    fprintf(fid, '%f', Sum4baseLI_value(sslp));fprintf(fid, '\t');
    fprintf(fid, '%f', Sum4baseLI_vn(sslp));fprintf(fid, '\t');
    fprintf(fid, '%d', Sum4VX_left(sslp));fprintf(fid, '\t');
    fprintf(fid, '%d', Sum4VX_right(sslp));fprintf(fid, '\t');
    fprintf(fid, '%f', Sum4value_left(sslp));fprintf(fid, '\t');
    fprintf(fid, '%f', Sum4value_right(sslp));fprintf(fid, '\t');
    fprintf(fid, '%f', Sum4AveVoxel_left(sslp));fprintf(fid, '\t');
    fprintf(fid, '%f', Sum4AveVoxel_right(sslp));fprintf(fid, '\n');

  end    % end of sslp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ************** Main Loop End **************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% write explanations
fprintf(fid, '\n');
fprintf(fid, '%s', 'Notes.');fprintf(fid, '\n');
fprintf(fid, '%s', 'AveLI: average lateralization index by voxel values with a weighting using all positive values (>0) (recommended)');fprintf(fid, '\n');
fprintf(fid, '%s', 'AveLI(voxel): index by voxel numbers computed using the AveLI algorithm');fprintf(fid, '\n');
fprintf(fid, '%s', 'baseLI: index by the simple summation (no weighting) of all positive voxel values (>0)');fprintf(fid, '\n');
fprintf(fid, '%s', 'baseLI(voxel): index by all voxel numbers with positive values (>0)');fprintf(fid, '\n');
fprintf(fid, '%s', 'Lt voxels: number of voxels with positive values (>0) in the left ROI');fprintf(fid, '\n');
fprintf(fid, '%s', 'Rt voxels: number of voxels with positive values (>0) in the right ROI');fprintf(fid, '\n');
fprintf(fid, '%s', 'Lt sum values: summation of all positive values (>0) in the left ROI');fprintf(fid, '\n');
fprintf(fid, '%s', 'Rt sum values: summation of all positive values (>0) in the right ROI');fprintf(fid, '\n');
fprintf(fid, '%s', 'Lt mean value: average of all positive values (>0) in the left ROI');fprintf(fid, '\n');
fprintf(fid, '%s', 'Rt mean value: average of all positive values (>0) in the right ROI');fprintf(fid, '\n');
fprintf(fid, '%s', '*NaN (if any) may mean that the ROI(s) do(es) not have voxels with positive values');fprintf(fid, '\n');


  fclose(fid);


% ****** Display results on the monitor ******

% Header of the log 
  disp(' ');disp(' ');
  disp('***** Results (details in the text file) *****');

% Information of specified files
  disp('>> Information of files you specified <<');
  disp('ROI files');
    x = ['Left ROI: ',leftmaskf];
  disp(x);
    x = ['Right ROI: ',rightmaskf];
  disp(x);
  disp('Data files');

  for i = 1:loopn_fff
      x = [num2str(i),': ',fff(i,:)];
    disp(x);
  end

  disp(' ');

    x = ['SS	AveLI	Lt voxels	Rt voxels'];
  disp(x);

  for i = 1:loopn_fff
      x = [num2str(i),'	',num2str(Sum4AveLI_value(i)),'	',num2str(Sum4VX_left(i)),'	',num2str(Sum4VX_right(i))];
    disp(x);
  end

  disp(' ');

  disp('Notes.');
  disp('AveLI    : average lateralization index by all positive voxel values (recommended)');
  disp('Lt voxels: number of voxels with positive values (>0) in the left ROI');
  disp('Rt voxels: number of voxels with positive values (>0) in the right ROI');
  disp('*NaN (if any) may mean that the ROI(s) do(es) not have voxels with positive values');


% ending display
  fprintf('\n')
  disp('   Good luck with your analysis!')
  fprintf('\n')

% save;

toc;

% END of script
