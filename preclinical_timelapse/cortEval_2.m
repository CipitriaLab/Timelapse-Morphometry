
% CORTICAL EVALUATION, Version 2
% Executes time-lapse morphometry evaluation of cortical bone. 
% Method is optimized for mouse tibiae

% Conceptualized by: Bettina M. Willie and Annette Birkhold.
% Written by: Annette Birkhold in 2011
% Functionalized by: Maximillian Rummler and Isabela Vitienes in 2019

% When used, please use the following reference:
% Birkhold A, Razi H, Duda GN, Weinkamer R, Checa S, Willie BM,
% Mineralizing surface is the main target of mechanical stimulation
% independent of age: 3D Dynamic in vivo Morphometry, Bone, 66:15-25, 2014

% LICENSE
% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation, either version 2 of the License, or (at your option) any later version.
% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
% For more details: https://www.gnu.org/licenses/gpl-3.0.html.

function [] = cortEval_2(inputfiles, paths, thresh, morph, days, is_short, is_mat, res_folder, reso)
tic
% INPUT
% inputfiles - nx1 cell array, ordered according to chronology of n scans
% paths - nx1 cell array, ordered as files
% thresh - cortical thresholds, nx1 array of doubles where n is no. of scans
% morph - 4x1 double array with morphometric radii
% is_short - boolean indicating if short report. if 0, long report
% is_mat - boolean indicating if output is mat file. if 0, tiff file
% res_folder - 2x1 cell array with paths for report (1) and output (2)
% reso - image resolution in um (i.e. this is the length of a single side
% of a voxel)

% Visualization of function progress
uf = uifigure;
name = uf.Name;
uf.Name = 'Running Evaluation';
uf.InnerPosition = [500 500 400 75];
progress = uiprogressdlg(uf,'Message','Please wait...','Indeterminate','on');

% Get files
reffile = inputfiles{1};
path = paths{1};
level = thresh(1);
files = inputfiles(2:end); % excludes ref file
path_2 = paths{2}; % excludes ref path

% create morphometric objects
R_1 = morph{1};
R_2 = morph{2};
R_3 = morph{3};
R_4 = morph{4};

disc_1 = strel('disk',R_1,8);
disc_2 = strel('disk',R_2,8);
disc_3 = strel('disk',R_3,8);
disc_4 = strel('disk',R_4,8);

a_1 = strrep(reffile, '.tif', '');
b_1 = isequal(a_1, reffile);
c_1 = strrep(a_1, '.mat', '');

a_2 = strrep(files, '.tif', '');
b_2 = isequal(a_2, reffile);
c_2 = strrep(a_2, '.mat', '');

% Open ref file
if b_1 == 1
    IMG_day0 = importdata([path reffile]);
    Image_day0(:,:,:) = IMG_day0(1,:,:,:);
else % input is tiff, slice-wise opening
    IMG_day0 = Tiff([path reffile], 'r');
    info = imfinfo([path reffile]);
    s = size(info);
    s = s(1);
    r1 = read(IMG_day0);
    S = size(r1);
    m = S(1);
    n = S(2);
    Image_day0 = zeros(m, n, s);
    for i = 1:s-1
        Image_day0(:,:,i) = read(IMG_day0);
        nextDirectory(IMG_day0);
    end
    Image_day0 = int16(Image_day0);
end

% Filtering ref file
Image_0 = smooth3(Image_day0,'gaussian',[3 3 3],0.65);
Image_tresh0 = Image_0 >= level;
clear IMG_day0 & Image_day0 & Image_0;

% Convert ref file into 8bit integer
Image_t0 = uint8(Image_tresh0);
clear Image_tresh0;

% File-wise evaluation
for m = 1:length(files) % non-baseline scans
    
    % Open file X
    if b_1 == 1
        IMG_dayX = importdata([path_2 char(files{m})]);
        Image_dayX(:,:,:) = IMG_dayX(1,:,:,:);
    else % input is tiff
        IMG_dayX = Tiff([path_2 char(files{m})], 'r');
        infox = imfinfo([path reffile]);
        sx = size(infox);
        sx = sx(1);
        rx = read(IMG_dayX);
        Sx = size(rx);
        nx = Sx(2);
        mx = Sx(1);
        Image_dayX = zeros(mx, nx, sx);
        for i = 1:sx-1
            Image_dayX(:,:,i) = read(IMG_dayX);
            nextDirectory(IMG_dayX);
        end
    end
    
    % Filtering/thresholding file X
    Image_X = smooth3(Image_dayX,'gaussian',[3 3 3],0.65);
    Image_treshX = Image_X >= thresh(m+1);
    clear IMG_dayX & Image_dayX & Image_X;
    
    % Convert file X to 8-bit 
    Image_tX = uint8(Image_treshX);
    clear Image_treshX;
    
    
    % ++++++++++++++++++++++++++++++++++++
    % +++++++CALCULATION OF VOLUMES+++++++
    % ++++++++++++++++++++++++++++++++++++
    
    % resorbed bone = volumes only in day 0
    resorbed_bone_X = (Image_tX == 0) & (Image_t0 == 1);
    resorbed_bone = uint8(resorbed_bone_X);
    
    % formed bone = volumes only in day X
    formed_bone_X = (Image_tX == 1) & (Image_t0 == 0);
    formed_bone = uint8(formed_bone_X);
    
    % constant bone = volumes in both datasets
    constant_bone_X = Image_t0 - resorbed_bone;
    constant_bone = uint8(constant_bone_X);
    
    % Creating labels
    Result_volumes_inclsurf_labels = constant_bone + formed_bone*2 + resorbed_bone*3;
    Result_volumes_inclsurf = constant_bone + formed_bone + resorbed_bone;
  
    
    % ++++++++++++++++++++++++++++++++++++
    % +++++++++++SEGMENTATION+++++++++++++
    % ++++++++++++++++++++++++++++++++++++
    
    % Perform slicewise segmentation on union of both data sets.
   
    [I_x,I_y,I_z] = size(Result_volumes_inclsurf);
    Endosteal_both = zeros(I_x,I_y,I_z);
    Periosteal_both = zeros(I_x,I_y,I_z);
    Cortical_both = zeros(I_x,I_y,I_z);
    
    for n = 1:I_z % slicewise segmentation
        I_bin = Result_volumes_inclsurf(:,:,n);
        % Dilation followed by closure to smoothen and close holes/lesions
        % propagating through cortex
        I2 = imdilate(I_bin,disc_1);
        I3 = imclose(I2,disc_2);
        % Filling cortical pores and medullary canal
        I4 = imfill(I3,8,'holes');
        I5 = imerode(I4,disc_3);
        % Create mask of filled cortex without medullary canal
        I6 = (I_bin == 0) & (I5 == 1);
        % Convert to uint8
        I7 = uint8(I6);
        % Fill remaining gaps
        I8 = imclose(I7,disc_4);
        % Erode again to remove any cortical bone from the cancellous bone
        I9 = imerode(I8,disc_1);
        % Intersection of initial binarized image and inverted cancellous
        % mask to obtain a mask of the cortex
        I10 = (I_bin == 1) & (I9 == 0);
        % Convert to uint8
        cort_bone = uint8(I10);
        
        % Separating endosteum and periosteum
        
        % dilate first the cortical bone mask
        cort_filled = imdilate(cort_bone, disc_4);
        % filling the medullary canal
        I11 = imfill(cort_filled,'holes');
        % erode again with same disc
        cort_filled2 = imerode(I11, disc_4);
        % use shrink operation to shrink mask
        I12 = bwmorph(cort_filled2,'shrink',3);
        % Endosteum
        I13 = (cort_bone == 1) & (I12 == 1);
        % Periosteum
        I14 = (cort_bone == 1) & (I12 == 0);
        I15 = uint8(I13);
        I16 = uint8(I14);
        Endosteal_both(:,:,n) = I15(:,:);
        Periosteal_both(:,:,n) = I16(:,:);
        Cortical_both(:,:,n) = cort_bone(:,:);
        clear n & I_bin & I2 & I3 & I4 & I5 & I6 & I7 & I8 & I9 & I10 & I11 & I12 & I13 & I14
        clear I15 & I16 & cort_bone & cort_filled & cort_filled2
    end
    clear I_x & I_y & I_z;
   
    
    % ++++++++++++++++++++++++++++++++++++
    % ++++++++++++MORPHOMETRY+++++++++++++
    % ++++++++++++++++++++++++++++++++++++
    
    % surface voxels:
    [I_x,I_y,I_z] = size(Cortical_both);
    Image_surf_labels = zeros(I_x,I_y,I_z);
    for n = 1:I_z %slicewise
        I_bin = Cortical_both(:,:,n);
        I2 = bwmorph(I_bin,'remove'); % deletes everything but surface voxels
        Image_surf_labels(:,:,n) = I2(:,:);
        clear I_bin & I2 & n;
    end
    clear I_x & I_y & I_z;
    
    % cortex without surface
    Image_without_surf_logic = (Cortical_both == 1) & (Image_surf_labels == 0);
    Image_without_surf = uint8(Image_without_surf_logic);
    [I_x,I_y,I_z] = size(Image_without_surf);
    Image_labels = zeros(I_x,I_y,I_z);
    for n = 1:I_z
        I_bin =  Image_without_surf(:,:,n);
        I2 = bwmorph(I_bin,'remove'); % deletes everything but surface voxels
        Image_labels(:,:,n) = I2(:,:); % 1st layer voxels
        clear I_bin & I2 & n;
    end
    clear I_x & I_y & I_z;
    
    % calculate SA of formed bone:
    % surface area of constant bone, grow 1, check if it´s newly formed
    [I_x,I_y,I_z] = size(Image_without_surf);
    Image_surf_const_labels = zeros(I_x,I_y,I_z);
    Constant_bone_volume1 = (Image_without_surf == 1) & (constant_bone == 1);
    Constant_bone_volume = uint8(Constant_bone_volume1);
    for n = 1:I_z
        I_bin = Constant_bone_volume(:,:,n);
        disc_5 = strel('disk',1,0);
        I2 = imdilate(I_bin,disc_5);
        I3 = bwmorph(I2,'remove');
        Image_surf_const_labels(:,:,n) = I3(:,:);
        clear I_bin & I2 & n & I3;
    end
    clear I_x & I_y & I_z;
    Image_surf_formed_labels1 =  (Image_surf_const_labels == 1) & (formed_bone == 1);
    Image_surf_formed_labels = uint8(Image_surf_formed_labels1);
    
    %+++++++++++++++++get surface of formation/resorption sites++++++++++++++++
    % whole bone
    const_bone_surf_voxels = (Image_labels == 1) & (constant_bone == 1);
    %     const_bone_surf_voxel = uint8(const_bone_surf_voxels);
    formed_bone_surf_voxels = Image_surf_formed_labels;
    resorbed_bone_surf_voxels = (Image_labels == 1) & (resorbed_bone == 1);
    % endosteal
    endostealArea = (Image_labels == 1) & (Endosteal_both == 1);
    endosteal_const_bone_surf_voxels = (endostealArea == 1) & (constant_bone == 1);
    endosteal_formed_bone_surf_voxels = (endostealArea == 1) & (formed_bone == 1);
    endosteal_resorbed_bone_surf_voxels = (endostealArea == 1) & (resorbed_bone == 1);
    %     endosteal_const_bone_surf_voxel = uint8(endosteal_const_bone_surf_voxels);
    % periosteal
    periostealArea = (Image_labels == 1) & (Periosteal_both == 1);
    periosteal_const_bone_surf_voxels = (periostealArea == 1) & (constant_bone == 1);
    periosteal_formed_bone_surf_voxels = (periostealArea == 1) & (formed_bone == 1);
    periosteal_resorbed_bone_surf_voxels = (periostealArea == 1) & (resorbed_bone == 1);
    %     periosteal_const_bone_surf_voxel = uint8(periosteal_const_bone_surf_voxels);
    
    % get amount of surface voxels
    % whole bone
    const_bone_surf_voxels_nr = nnz(const_bone_surf_voxels);
    formed_bone_surf_voxels_nr = nnz(formed_bone_surf_voxels);
    resorbed_bone_surf_voxels_nr = nnz(resorbed_bone_surf_voxels);
    % endosteal
    endosteal_const_bone_surf_voxels_nr = nnz (endosteal_const_bone_surf_voxels);
    endosteal_formed_bone_surf_voxels_nr = nnz (endosteal_formed_bone_surf_voxels);
    endosteal_resorbed_bone_surf_voxels_nr = nnz (endosteal_resorbed_bone_surf_voxels);
    % periosteal
    periosteal_const_bone_surf_voxels_nr = nnz (periosteal_const_bone_surf_voxels);
    periosteal_formed_bone_surf_voxels_nr = nnz (periosteal_formed_bone_surf_voxels);
    periosteal_resorbed_bone_surf_voxels_nr = nnz (periosteal_resorbed_bone_surf_voxels);
    
    % convert voxel into area
    rArea = reso^2;
    % whole bone
    const_bone_surf_area = const_bone_surf_voxels_nr * rArea;
    formed_bone_surf_area = formed_bone_surf_voxels_nr * rArea;
    resorbed_bone_surf_area = resorbed_bone_surf_voxels_nr * rArea;
    % endosteal
    endosteal_const_bone_surf_area = endosteal_const_bone_surf_voxels_nr * rArea;
    endosteal_formed_bone_surf_area = endosteal_formed_bone_surf_voxels_nr * rArea;
    endosteal_resorbed_bone_surf_area = endosteal_resorbed_bone_surf_voxels_nr * rArea;
    % periosteal
    periosteal_const_bone_surf_area = periosteal_const_bone_surf_voxels_nr * rArea;
    periosteal_formed_bone_surf_area = periosteal_formed_bone_surf_voxels_nr * rArea;
    periosteal_resorbed_bone_surf_area = periosteal_resorbed_bone_surf_voxels_nr * rArea;
    
    % total area day 0
    % whole bone
    % edits: added formed bone surface area to total area calculation
    % based on LAIA's correction
    totalArea = const_bone_surf_area + resorbed_bone_surf_area + formed_bone_surf_area;
    const_SATA = const_bone_surf_area/totalArea;
    formed_SATA = formed_bone_surf_area/totalArea;
    resorbed_SATA = resorbed_bone_surf_area/totalArea;
    clear totalArea;
    % endosteal
    endosteal_totalArea = endosteal_const_bone_surf_area + endosteal_resorbed_bone_surf_area + endosteal_formed_bone_surf_area;
    endosteal_const_SATA = endosteal_const_bone_surf_area/endosteal_totalArea;
    endosteal_formed_SATA = endosteal_formed_bone_surf_area/endosteal_totalArea;
    endosteal_resorbed_SATA = endosteal_resorbed_bone_surf_area/endosteal_totalArea;
    clear endosteal_totalArea
    % periosteal
    periosteal_totalArea = periosteal_const_bone_surf_area +  periosteal_resorbed_bone_surf_area + periosteal_formed_bone_surf_area;
    periosteal_const_SATA = periosteal_const_bone_surf_area/periosteal_totalArea;
    periosteal_formed_SATA = periosteal_formed_bone_surf_area/periosteal_totalArea;
    periosteal_resorbed_SATA = periosteal_resorbed_bone_surf_area/periosteal_totalArea;
    clear endosteal_totalArea
    
    % calculate volumes
    % whole bone
    constant_bone_voxels_all_logic = (Image_without_surf==1) & (constant_bone==1);
    formed_bone_voxels_all_logic = (Image_without_surf==1) & (formed_bone==1);
    resorbed_bone_voxels_all_logic = (Image_without_surf==1) & (resorbed_bone==1);
    % endosteal
    endostealVolume = (Image_without_surf == 1) & (Endosteal_both == 1);
    %     endosteal_constant_bone_voxels_all_logic = (endostealVolume==1) & (constant_bone==1);
    endosteal_formed_bone_voxels_all_logic = (endostealVolume==1) & (formed_bone==1);
    endosteal_resorbed_bone_voxels_all_logic = (endostealVolume==1) & (resorbed_bone==1);
    % periosteal
    periostealVolume = (Image_labels == 1) & (Periosteal_both == 1);
    %     periosteal_constant_bone_voxels_all_logic = (periostealVolume==1) & (constant_bone==1);
    periosteal_formed_bone_voxels_all_logic = (periostealVolume==1) & (formed_bone==1);
    periosteal_resorbed_bone_voxels_all_logic = (periostealVolume==1) & (resorbed_bone==1);
    
    % number of voxel in the total volume without the surface
    % whole bone
    formed_bone_voxels_all = nnz(formed_bone_voxels_all_logic);
    constant_bone_voxels_all = nnz(constant_bone_voxels_all_logic);
    resorbed_bone_voxels_all = nnz(resorbed_bone_voxels_all_logic);
    clear resorbed_bone_voxels_all_logic & constant_bone_voxels_all_logic & formed_bone_voxels_all_logic;
    % endosteal
    endosteal_formed_bone_voxels_all = nnz(endosteal_formed_bone_voxels_all_logic);
    %     endosteal_constant_bone_voxels_all = nnz(endosteal_constant_bone_voxels_all_logic);
    endosteal_resorbed_bone_voxels_all = nnz(endosteal_resorbed_bone_voxels_all_logic);
    clear endosteal_resorbed_bone_voxels_logic & endosteal_constant_bone_voxels_logic & endosteal_formed_bone_voxels_logic
    % periosteal
    periosteal_formed_bone_voxels_all = nnz(periosteal_formed_bone_voxels_all_logic);
    %     periosteal_constant_bone_voxels_all = nnz(periosteal_constant_bone_voxels_all_logic);
    periosteal_resorbed_bone_voxels_all = nnz(periosteal_resorbed_bone_voxels_all_logic);
    clear periosteal_resorbed_bone_voxels_logic & periosteal_constant_bone_voxels_logic & periosteal_formed_bone_voxels_logic
    
    % number of voxel in the total volume with the surface -> half of
    % the surface voxel added due to partial volume effect -> the
    % original surface is not added!
    % whole bone
    constant_bone_voxels = constant_bone_voxels_all + const_bone_surf_voxels_nr*0.5;
    formed_bone_voxels = formed_bone_voxels_all + formed_bone_surf_voxels_nr*0.5;
    resorbed_bone_voxels = resorbed_bone_voxels_all + resorbed_bone_surf_voxels_nr*0.5;
    clear resorbed_bone_voxels_all & constant_bone_voxels_all & formed_bone_voxels_all;
    % endosteal
    %     endosteal_constant_bone_voxels = endosteal_constant_bone_voxels_all + endosteal_const_bone_surf_voxels_nr*0.5;
    endosteal_formed_bone_voxels = endosteal_formed_bone_voxels_all + endosteal_formed_bone_surf_voxels_nr*0.5;
    endosteal_resorbed_bone_voxels = endosteal_resorbed_bone_voxels_all + endosteal_resorbed_bone_surf_voxels_nr*0.5;
    clear endosteal_resorbed_bone_voxels_all & endosteal_constant_bone_voxels_all & endosteal_formed_bone_voxels_all
    % periosteal
    %     periosteal_constant_bone_voxels = periosteal_constant_bone_voxels_all + periosteal_const_bone_surf_voxels_nr*0.5;
    periosteal_formed_bone_voxels = periosteal_formed_bone_voxels_all + periosteal_formed_bone_surf_voxels_nr*0.5;
    periosteal_resorbed_bone_voxels = periosteal_resorbed_bone_voxels_all + periosteal_resorbed_bone_surf_voxels_nr*0.5;
    clear periosteal_resorbed_bone_voxels_all & periosteal_constant_bone_voxels_all & periosteal_formed_bone_voxels_all
    
    % convert voxel into volume
    rVol = reso^3;
    % whole bone
    constant_bone_volume = constant_bone_voxels * rVol;
    formed_bone_volume = formed_bone_voxels * rVol;
    resorbed_bone_volume = resorbed_bone_voxels * rVol;
    clear resorbed_bone_voxels & formed_bone_voxels & constant_bone_voxels;
    % endosteal
    %     endosteal_constant_bone_volume = endosteal_constant_bone_voxels * rVol;
    endosteal_formed_bone_volume = endosteal_formed_bone_voxels * rVol;
    endosteal_resorbed_bone_volume = endosteal_resorbed_bone_voxels * rVol;
    clear endosteal_resorbed_bone_voxels & endosteal_formed_bone_voxels & endosteal_constant_bone_voxels
    % periosteal
    %     periosteal_constant_bone_volume = periosteal_constant_bone_voxels * rVol;
    periosteal_formed_bone_volume = periosteal_formed_bone_voxels * rVol;
    periosteal_resorbed_bone_volume = periosteal_resorbed_bone_voxels * rVol;
    clear periosteal_resorbed_bone_voxels & periosteal_formed_bone_voxels & periosteal_constant_bone_voxels
    
    %VOLUME/total volume % total bone volume day 0
    % whole bone
    totalVolume     = constant_bone_volume + resorbed_bone_volume; % true volume at reference day!
    constant_BVTV   = constant_bone_volume / totalVolume;
    formed_BVTV  	= formed_bone_volume / totalVolume;
    resorbed_BVTV   = resorbed_bone_volume / totalVolume;
    % endosteal
    endosteal_formed_BVTV = endosteal_formed_bone_volume / totalVolume;
    endosteal_resorbed_BVTV = endosteal_resorbed_bone_volume / totalVolume;
    % periosteal
    periosteal_formed_BVTV = periosteal_formed_bone_volume / totalVolume;
    periosteal_resorbed_BVTV = periosteal_resorbed_bone_volume / totalVolume;
    clear totalVolume;
    
    % average height of volumes
    % whole bone
    formed_height = formed_bone_volume / formed_bone_surf_area;
    resorb_height = resorbed_bone_volume / resorbed_bone_surf_area;
    % endosteal
    endosteal_formed_height = endosteal_formed_bone_volume / endosteal_formed_bone_surf_area;
    endosteal_resorb_height = endosteal_resorbed_bone_volume / endosteal_resorbed_bone_surf_area;
    % periosteal
    periosteal_formed_height = periosteal_formed_bone_volume / periosteal_formed_bone_surf_area;
    periosteal_resorb_height = periosteal_resorbed_bone_volume / periosteal_resorbed_bone_surf_area;
    
    
    % following the average thickness of formed and resorbed bone is calculated
    % mean thickness formed bone
    % whole bone
    I_formed_without_surf = double(Image_without_surf).* double(formed_bone);
    I_inv = -(I_formed_without_surf - 1);
    Dist = bwdist(I_inv,'euclidean'); % euclidian distance map of image
    Ult_E = bwulterode(Dist,'euclidean',6); % ultimate erosion of the distance map -> gets the local maxima -> point where most bone was formed
    E = double(Ult_E);
    dist_formed = 2 * double(Dist) - 1;
    loc_max_formed = E .* dist_formed .* reso; % Voxelsize
    % local maxima*value at local maxima = E*Dist; µm
    mean_formed =  mean (loc_max_formed(loc_max_formed~=0));
    std_formed = std (loc_max_formed(loc_max_formed~=0));
    max_formed = max(loc_max_formed(loc_max_formed~=0));
    clear I_formed_without_surf & I_inv & Dist & Ult_E & E;
    % endosteal
    I_formed_without_surf = double(Image_without_surf).* double(formed_bone) & (Endosteal_both == 1);
    I_inv = -(I_formed_without_surf - 1);
    Dist = bwdist(I_inv,'euclidean'); % euclidian distance map of image
    Ult_E = bwulterode(Dist,'euclidean',6); % ultimate erosion of the distance map -> gets the local maxima -> point where most bone was formed
    E = double(Ult_E);
    dist_formedEndosteal = 2 * double(Dist) - 1;
    loc_max_formedEndosteal = E .* dist_formedEndosteal .* reso;
    % local maxima*value at local maxima = E*Dist; µm
    mean_formedEndosteal = mean (loc_max_formedEndosteal(loc_max_formedEndosteal~=0));
    std_formedEndosteal = std (loc_max_formedEndosteal(loc_max_formedEndosteal~=0));
    max_formedEndosteal = max(loc_max_formedEndosteal(loc_max_formedEndosteal~=0));
    clear I_formed_without_surf & I_inv & Dist & Ult_E & E;
    % periosteal
    I_formed_without_surf = double(Image_without_surf).* double(formed_bone) & (Periosteal_both == 1);
    I_inv = -(I_formed_without_surf - 1);
    Dist = bwdist(I_inv,'euclidean'); % euclidian distance map of image
    Ult_E = bwulterode(Dist,'euclidean',6); % ultimate erosion of the distance map -> gets the local maxima -> point where most bone was formed
    E = double(Ult_E);
    dist_formedPeriosteal = 2 * double(Dist) - 1;
    loc_max_formedPeriosteal = E .* dist_formedPeriosteal .* reso;
    % local maxima*value at local maxima = E*Dist; µm
    mean_formedPeriosteal = mean (loc_max_formedPeriosteal(loc_max_formedPeriosteal~=0));
    std_formedPeriosteal = std (loc_max_formedPeriosteal(loc_max_formedPeriosteal~=0));
    max_formedPeriosteal = max(loc_max_formedPeriosteal(loc_max_formedPeriosteal~=0));
    clear I_formed_without_surf & I_inv & Dist & Ult_E & E;
    
    % mean thickness resorbed bone
    % whole bone
    I_resorbed_without_surf = double(Image_without_surf).* double(resorbed_bone);
    I_inv = -(I_resorbed_without_surf - 1);
    Dist = bwdist(I_inv,'euclidean'); % euclidian distance map of image
    Ult_E = bwulterode(Dist,'euclidean',6); % ultimate erosion of the distance map -> gets the local maxima -> point where most bone was formed
    E = double(Ult_E);
    dist_resorbed = 2 * double(Dist) - 1;
    loc_max_resorbed = E .* dist_resorbed .* reso;
    % local maxima*value at local maxima = E*Dist; µm
    mean_resorbed =  mean (loc_max_resorbed(loc_max_resorbed~=0));
    std_resorbed = std (loc_max_resorbed(loc_max_resorbed~=0));
    max_resorbed = max(loc_max_resorbed(loc_max_resorbed~=0));
    clear I_resorbed_without_surf & I_inv & Dist & Ult_E & E;
    % endosteal bone
    I_resorbed_without_surf = double(Image_without_surf).* double(resorbed_bone) & (Endosteal_both == 1);
    I_inv = -(I_resorbed_without_surf - 1);
    Dist = bwdist(I_inv,'euclidean'); % euclidian distance map of image
    Ult_E = bwulterode(Dist,'euclidean',6); % ultimate erosion of the distance map -> gets the local maxima -> point where most bone was formed
    E = double(Ult_E);
    dist_resorbedEndosteal = 2 * double(Dist) - 1;
    loc_max_resorbedEndosteal = E .* dist_resorbedEndosteal .* reso;
    % local maxima*value at local maxima = E*Dist; µm
    mean_resorbedEndosteal = mean (loc_max_resorbedEndosteal(loc_max_resorbedEndosteal~=0));
    std_resorbedEndosteal = std (loc_max_resorbedEndosteal(loc_max_resorbedEndosteal~=0));
    max_resorbedEndosteal = max(loc_max_resorbedEndosteal(loc_max_resorbedEndosteal~=0));
    clear I_resorbed_without_surf & I_inv & Dist & Ult_E & E;
    % periosteal bone
    I_resorbed_without_surf = double(Image_without_surf).* double(resorbed_bone) & (Periosteal_both == 1);
    I_inv = -(I_resorbed_without_surf - 1);
    Dist = bwdist(I_inv,'euclidean'); % euclidian distance map of image
    Ult_E = bwulterode(Dist,'euclidean',6); % ultimate erosion of the distance map -> gets the local maxima -> point where most bone was formed
    E = double(Ult_E);
    dist_resorbedPeriosteal = 2 * double(Dist) - 1;
    loc_max_resorbedPeriosteal = E .* dist_resorbedPeriosteal .* reso;
    % local maxima*value at local maxima = E*Dist; µm
    mean_resorbedPeriosteal = mean (loc_max_resorbedPeriosteal(loc_max_resorbedPeriosteal~=0));
    std_resorbedPeriosteal = std (loc_max_resorbedPeriosteal(loc_max_resorbedPeriosteal~=0));
    max_resorbedPeriosteal = max(loc_max_resorbedPeriosteal(loc_max_resorbedPeriosteal~=0));
    clear I_resorbed_without_surf & I_inv & Dist & Ult_E & E;
    
    
    % 3D - MAR 
    day = cell2mat(days(m+1));
    % whole bone
    MAR_3D = mean_formed/day;
    MRR_3D = mean_resorbed/day;
    % endosteal
    endosteal_MAR_3D = mean_formedEndosteal / day;
    endosteal_MRR_3D = mean_resorbedEndosteal / day;
    % periosteal
    periosteal_MAR_3D = mean_formedPeriosteal / day;
    periosteal_MRR_3D = mean_resorbedPeriosteal / day;    
    
    % Volume of bone formed/resorbed per time unit
    % whole bone
    BFR_3D = formed_bone_volume/day; % bone volume formed per day
    BRR_3D = resorbed_bone_volume/day;
    % endosteal
    endosteal_BFR_3D = endosteal_formed_bone_volume/day; % bone volume formed per day
    endosteal_BRR_3D = endosteal_resorbed_bone_volume/day;
    % periosteal
    periosteal_BFR_3D = periosteal_formed_bone_volume/day; % bone volume formed per day
    periosteal_BRR_3D = periosteal_resorbed_bone_volume/day;
    
    % Write volumetric results
    filename = char(c_2(m));
    if is_short == 1
        Parameters = {'Threshold'; 'R_1'; 'R_2'; 'R_3'; 'R_4'; 'Resolution'; 'Ec.MV/BV'; 'Ec.EV/BV'; 'Ec.MS/BS'; 'Ec.ES/BS'; 'Ec.MAR3D'; 'Ec.MRR3D'; 'Ec.BFR3D'; 'Ec.BRR3D'; 'Ec.MTh'; 'Ec.ED'; 'Ps.MV/BV'; 'Ps.EV/BV'; 'Ps.MS/BS'; 'Ps.ES/BS'; 'Ps.MAR3D'; 'Ps.MRR3D'; 'Ps.BFR3D'; 'Ps.BRR3D'; 'Ps.MTh'; 'Ps.ED'; 'T.MV/BV'; 'T.EV/BV'; 'T.MS/BS'; 'T.ES/BS'; 'T.MAR3D'; 'T.MRR3D'; 'T.BFR3D'; 'T.BRR3D'; 'T.MTh'; 'T.ED'};
        Values = [thresh(m+1); R_1; R_2; R_3; R_4; reso; endosteal_formed_BVTV; endosteal_resorbed_BVTV; endosteal_formed_SATA; endosteal_resorbed_SATA; endosteal_MAR_3D; endosteal_MRR_3D; endosteal_BFR_3D; endosteal_BRR_3D; mean_formedEndosteal; mean_resorbedEndosteal; periosteal_formed_BVTV; periosteal_resorbed_BVTV; periosteal_formed_SATA; periosteal_resorbed_SATA; periosteal_MAR_3D; periosteal_MRR_3D; periosteal_BFR_3D; periosteal_BRR_3D; mean_formedPeriosteal; mean_resorbedPeriosteal; formed_BVTV; resorbed_BVTV; formed_SATA; resorbed_SATA; MAR_3D; MRR_3D; BFR_3D; BRR_3D; formed_height; resorb_height];        
        resultsfile = strcat(char(res_folder(1)), '\', filename, '_cort_VolResults_short.xlsx');
    else
        
        Parameters = {'Threshold'; 'R_1'; 'R_2'; 'R_3'; 'R_4'; 'Resolution'; 'Ec.MV'; 'Ec.EV' ; 'Ec.MV/BV'; 'Ec.EV/BV'; 'Ec.CS'; 'Ec.MS'; 'Ec.ES'; 'Ec.CS/BS'; 'Ec.MS/BS'; 'Ec.ES/BS'; 'Ec.MAR3D'; 'Ec.MRR3D'; 'Ec.BFR3D'; 'Ec.BRR3D'; 'Ec.MTh'; 'Ec.MTh.sd'; 'Ec.MTh.max'; 'Ec.ED'; 'Ec.ED.sd'; 'Ec.ED.max'; 'Ec.Md'; 'Ec.Ed';'Pt.MV'; 'Ps.EV' ; 'Ps.MV/BV'; 'Ps.EV/BV'; 'Ps.CS'; 'Ps.MS'; 'Ps.ES'; 'Ps.CS/BS'; 'Ps.MS/BS'; 'Ps.ES/BS'; 'Ps.MAR3D'; 'Ps.MRR3D'; 'Ps.BFR3D'; 'Ps.BRR3D'; 'Ps.MTh'; 'Ps.MTh.sd'; 'Ps.MTh.max'; 'Ps.ED'; 'Ps.ED.sd'; 'Ps.ED.max'; 'Ps.Md'; 'Ps.Ed'; 'T.CV'; 'T.MV'; 'T.EV'; 'T.CV/BV'; 'T.MV/BV'; 'T.EV/BV'; 'T.CS'; 'T.MS'; 'T.ES'; 'T.CS/BS'; 'T.MS/BS'; 'T.ESBS'; 'T.MAR3D'; 'T.MRR3D'; 'T.BFR3D'; 'T.BRR3D'; 'T.MTh'; 'T.MTh.sd'; 'T.MTh.max'; 'T.ED'; 'T.ED.sd'; 'T.ED.max'; 'T.Md'; 'T.Ed'};
        Values = [thresh(m+1); R_1; R_2; R_3; R_4; reso; endosteal_formed_bone_volume; endosteal_resorbed_bone_volume; endosteal_formed_BVTV; endosteal_resorbed_BVTV; endosteal_const_bone_surf_area; endosteal_formed_bone_surf_area; endosteal_resorbed_bone_surf_area; endosteal_const_SATA; endosteal_formed_SATA; endosteal_resorbed_SATA; endosteal_MAR_3D; endosteal_MRR_3D; endosteal_BFR_3D; endosteal_BRR_3D; mean_formedEndosteal; std_formedEndosteal; max_formedEndosteal; mean_resorbedEndosteal; std_resorbedEndosteal; max_resorbedEndosteal; endosteal_formed_height; endosteal_resorb_height; periosteal_formed_bone_volume; periosteal_resorbed_bone_volume; periosteal_formed_BVTV; periosteal_resorbed_BVTV; periosteal_const_bone_surf_area; periosteal_formed_bone_surf_area; periosteal_resorbed_bone_surf_area; periosteal_const_SATA; periosteal_formed_SATA; periosteal_resorbed_SATA; periosteal_MAR_3D; periosteal_MRR_3D; periosteal_BFR_3D; periosteal_BRR_3D; periosteal_formed_height; periosteal_resorb_height; mean_formedPeriosteal; std_formedPeriosteal; max_formedPeriosteal; mean_resorbedPeriosteal; std_resorbedPeriosteal; max_resorbedPeriosteal; constant_bone_volume; formed_bone_volume; resorbed_bone_volume; constant_BVTV; formed_BVTV; resorbed_BVTV; const_bone_surf_area; formed_bone_surf_area; resorbed_bone_surf_area; const_SATA; formed_SATA; resorbed_SATA; MAR_3D; MRR_3D; BFR_3D; BRR_3D; mean_formed; std_formed; max_formed; mean_resorbed; std_resorbed; max_resorbed; formed_height; resorb_height];       
        resultsfile = strcat(char(res_folder(1)), '', filename, '_cort_VolResults_long.xlsx');
    end
    T = table(Parameters, Values);
    writetable(T, resultsfile);
    
    
    % Vis
    if is_short == 1
        if is_mat == 1
            file_cortical_complete_VISres = strcat(filename, '_cortical_complete_ImageResults_short');
            path_cortical_complete_VISres = strcat(char(res_folder(2)), '/');
            savefile2 = ([path_cortical_complete_VISres file_cortical_complete_VISres]);
            save(savefile2, 'Result_volumes_inclsurf_labels' , 'Result_volumes_inclsurf');
        else % is tif
            Day_cell = num2str(day);
            path_cortical_complete_VISres = strcat(char(res_folder(2)), '/');
            
            RVISL = strcat(path_cortical_complete_VISres, 'Result_volumes_inclsurf_labels_', Day_cell, '.tif');      
            no_slices = size(Result_volumes_inclsurf_labels); no_slices = no_slices(3);
            for i = 1:no_slices
                imwrite(Result_volumes_inclsurf_labels(:,:,i), RVISL, 'WriteMode', 'append', 'Compression', 'none');
            end
            
            RVIS = strcat(path_cortical_complete_VISres, 'Result_volumes_inclsurf_', Day_cell, '.tif');
            no_slices = size(Result_volumes_inclsurf); no_slices = no_slices(3);
            for i = 1:no_slices
                imwrite(Result_volumes_inclsurf(:,:,i), RVIS, 'WriteMode', 'append', 'Compression', 'none');
            end
            
        end
    else
        if is_mat == 1
            file_cortical_complete_VISres = strcat(filename, '_cortical_complete_ImageResults_long');
            path_cortical_complete_VISres = strcat(char(res_folder(2)), '/');
            
            savefile2 = ([path_cortical_complete_VISres file_cortical_complete_VISres]);
            save(savefile2, 'Cortical_both', 'Image_t0', 'Image_tX', 'Result_volumes_inclsurf_labels' , 'Result_volumes_inclsurf');
            save(savefile2,  'Image_without_surf','-append' )
            save(savefile2,'constant_bone','formed_bone','resorbed_bone','-append' );
            save(savefile2, 'Endosteal_both' ,'Periosteal_both','-append')
            save(savefile2,'dist_formed','loc_max_formed','dist_resorbed','loc_max_resorbed','-append' );
        else % is tif
            Day_cell = num2str(day);
            path_cortical_complete_VISres = strcat(char(res_folder(2)),'\');
            
            corticalboth = strcat(path_cortical_complete_VISres, 'corticalBoth', Day_cell, '.tif');
            no_slices = size(Cortical_both); no_slices = no_slices(3);
            for i = 1:no_slices
                imwrite(Cortical_both(:,:,i), corticalboth, 'WriteMode', 'append', 'Compression', 'none');
            end
            
            it0 = strcat(path_cortical_complete_VISres, 'Image_t0_', Day_cell, '.tif');
            no_slices = size(Image_t0); no_slices = no_slices(3);
            for i = 1:no_slices
                imwrite(Image_t0(:,:,i), it0, 'WriteMode', 'append', 'Compression', 'none');
            end
            
            itX = strcat(path_cortical_complete_VISres, 'Image_tX_', Day_cell, '.tif');
            no_slices = size(Image_tX); no_slices = no_slices(3);
            for i = 1:no_slices
                imwrite(Image_tX(:,:,i), itX, 'WriteMode', 'append', 'Compression', 'none');
            end
            
            RVISL = strcat(path_cortical_complete_VISres, 'Result_volumes_inclsurf_labels_', Day_cell, '.tif');
            no_slices = size(Result_volumes_inclsurf_labels); no_slices = no_slices(3);
            for i = 1:no_slices
                imwrite(Result_volumes_inclsurf_labels(:,:,i), RVISL, 'WriteMode', 'append', 'Compression', 'none');
            end
            
            RVIS = strcat(path_cortical_complete_VISres, 'Result_volumes_inclsurf_', Day_cell, '.tif');
            no_slices = size(Result_volumes_inclsurf); no_slices = no_slices(3);
            for i = 1:no_slices
                imwrite(Result_volumes_inclsurf(:,:,i), RVIS, 'WriteMode', 'append', 'Compression', 'none');
            end
            
            imagewithoutsurf = strcat(path_cortical_complete_VISres, 'Image_without_surf_', Day_cell, '.tif');
            no_slices = size(Image_without_surf); no_slices = no_slices(3);
            for i = 1:no_slices
                imwrite(Image_without_surf(:,:,i), imagewithoutsurf, 'WriteMode', 'append', 'Compression', 'none');
            end
            
            constantbone = strcat(path_cortical_complete_VISres, 'constant_bone_', Day_cell, '.tif');
            no_slices = size(constant_bone); no_slices = no_slices(3);
            for i = 1:no_slices
                imwrite(constant_bone(:,:,i), constantbone, 'WriteMode', 'append', 'Compression', 'none');
            end
            
            formedbone = strcat(path_cortical_complete_VISres, 'formed_bone_', Day_cell, '.tif');
            no_slices = size(formed_bone); no_slices = no_slices(3);
            for i = 1:no_slices
                imwrite(formed_bone(:,:,i), formedbone, 'WriteMode', 'append', 'Compression', 'none');
            end
            
            resorbedbone = strcat(path_cortical_complete_VISres, 'resorbed_bone_', Day_cell, '.tif');
            no_slices = size(resorbed_bone); no_slices = no_slices(3);
            for i = 1:no_slices
                imwrite(resorbed_bone(:,:,i), resorbedbone, 'WriteMode', 'append', 'Compression', 'none');
            end
            
            Endostealboth = strcat(path_cortical_complete_VISres, 'Endosteal_both_', Day_cell, '.tif');
            no_slices = size(Endosteal_both); no_slices = no_slices(3);
            for i = 1:no_slices
                imwrite(Endosteal_both(:,:,i), Endostealboth, 'WriteMode', 'append', 'Compression', 'none');
            end
            
            Periostealboth = strcat(path_cortical_complete_VISres, 'Periosteal_both_', Day_cell, '.tif');
            no_slices = size(Periosteal_both); no_slices = no_slices(3);
            for i = 1:no_slices
                imwrite(Periosteal_both(:,:,i), Periostealboth, 'WriteMode', 'append', 'Compression', 'none');
            end
            
            distformed = strcat(path_cortical_complete_VISres, 'dist_formed_', Day_cell, '.tif');
            no_slices = size(dist_formed); no_slices = no_slices(3);
            for i = 1:no_slices
                imwrite(dist_formed(:,:,i), distformed, 'WriteMode', 'append', 'Compression', 'none');
            end
            
            locmaxformed = strcat(path_cortical_complete_VISres, 'loc_max_formed_', Day_cell, '.tif');
            no_slices = size(loc_max_formed); no_slices = no_slices(3);
            for i = 1:no_slices
                imwrite(loc_max_formed(:,:,i), locmaxformed, 'WriteMode', 'append', 'Compression', 'none');
            end
            
            distres = strcat(path_cortical_complete_VISres, 'dist_resorbed_', Day_cell, '.tif');
            no_slices = size(dist_resorbed); no_slices = no_slices(3);
            for i = 1:no_slices
                imwrite(dist_resorbed(:,:,i), distres, 'WriteMode', 'append', 'Compression', 'none');
            end
            
            locmaxres = strcat(path_cortical_complete_VISres, 'loc_max_resorbed_', Day_cell, '.tif');
            no_slices = size(loc_max_resorbed); no_slices = no_slices(3);
            for i = 1:no_slices
                imwrite(loc_max_resorbed(:,:,i), locmaxres, 'WriteMode', 'append', 'Compression', 'none');
            end
            
        end
    end

    
end
close(progress);
close(uf);
toc
end
