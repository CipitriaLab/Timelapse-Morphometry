function [] = Both_Eval(inputfiles, paths, thresh, morph, days, is_short, is_mat, res_folder, voxelsize, bone_main_axis,SegmentationSwitch, VisualizationSwitch,which_eval)
tic
% ======================================================================= %
%% 2022-03-10 Update
% Cleaning the code and implementing the new part from Sarah
% ======================================================================= %
% INPUT
% inputfiles - nx1 cell array, ordered according to chronology of n scans
% paths - nx1 cell array, ordered as files
% thresh - cortical threshold, double
% morph - 4x1 double array with morphometric radii
% is_short - boolean indicating if short report. if 0, long report
% is_mat - boolean indicating if output is mat file. if 0, tiff file
% res_folder - 2x1 cell array with paths for report (1) and output (2)

% Visualization of function progress
uf = uifigure;
name = uf.Name;
uf.Name = 'Running Evaluation';
uf.InnerPosition = [500 500 400 75];
progress = uiprogressdlg(uf,'Message','Please wait...','Indeterminate','on');

% Get files
reffile = inputfiles{1};
path = paths{1};
%% 2022-03-10 Update
% Having two separate threshold
level_0 = thresh(1,1);
level_x = thresh(2,1);
%%
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
%% 2022-03-10 Update
% switching x,y,z directions to y,z,x -> different position for
% segmentation
Image_day0 = switching_direction(Image_day0, bone_main_axis);

%% Filtering ref file
Image_0 = smooth3(Image_day0,'gaussian',[3 3 3],0.65);

%% 2022-03-10 Update
% Change Hub and Sarah
% Image_tresh0 = Image_0 >= level; old
Image_tresh0 = Image_0 >= level_0; % Use the new separate threshold
%%
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
    %% 2022-03-10 Update
    % Same as before, swaping the direction of the images.
    Image_dayX = switching_direction(Image_dayX, bone_main_axis);
    %% Filtering/thresholding file X
    Image_X = smooth3(Image_dayX,'gaussian',[3 3 3],0.65);
    %% 2022-03-10 Update
    % Change Hub and Sarah to use the two separate threshold.
    Image_treshX = Image_X >= level_x; %
    %%
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
    
    %% 2022-03-10 Update
    % Mod V1.3 Only having the constant and resorbed for the surface
    % evaluation of the surface of the reference bone at D0.
    Result_volumes_inclsurf_d0 = constant_bone + resorbed_bone;
    %%
    % ++++++++++++++++++++++++++++++++++++
    % +++++++++++SEGMENTATION+++++++++++++
    % ++++++++++++++++++++++++++++++++++++
    
    % Perform slicewise segmentation on union of both data sets.
    
    [I_x,I_y,I_z] = size(Result_volumes_inclsurf);
    switch which_eval
        case 'cort'
            Endosteal_both = zeros(I_x,I_y,I_z);
            Periosteal_both = zeros(I_x,I_y,I_z);
        case 'trab'
    end
    Bone_both = zeros(I_x,I_y,I_z); % can be for corticol or trab
    
    for n = 1:I_z % slicewise segmentation
        %% 2022-03-10 Update
        % In the app there is a switch to perform the segmentation or not.
        if strcmp(SegmentationSwitch,'On')
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
            switch which_eval
                case 'cort'
                    cort_bone = uint8(I10);
                case 'trab'
                    I11 = (I_bin == 1) & (I10 == 0);
                    % Convert to uint8
                    trab_bone = uint8(I11);
                    Bone_both(:,:,n) = uint8(trab_bone(:,:));
            end
        else
            % Do not do the segmentation.
            Bone_both = uint8(Result_volumes_inclsurf(:,:,n));
        end
        switch which_eval
            case 'cort'
                %% ONLY CORTICaL
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
                
                Bone_both(:,:,n) = cort_bone(:,:);
            case 'trab'
                % Do nothing here
        end
    end
    clear I_x & I_y & I_z;
    
    
    % ++++++++++++++++++++++++++++++++++++
    % ++++++++++++MORPHOMETRY+++++++++++++
    % ++++++++++++++++++++++++++++++++++++
    
    % surface voxels:
    [I_x,I_y,I_z] = size(Bone_both);
    Image_surf_labels = zeros(I_x,I_y,I_z);
    for n = 1:I_z %slicewise
        I_bin = Bone_both(:,:,n);
        I2 = bwmorph(I_bin,'remove'); % deletes everything but surface voxels
        Image_surf_labels(:,:,n) = I2(:,:);
        %         clear I_bin & I2 & n;
    end
    clear I_x & I_y & I_z;
    
    % cortex without surface
    Image_without_surf_logic = (Bone_both == 1) & (Image_surf_labels == 0);
    Image_without_surf = uint8(Image_without_surf_logic);
    [I_x,I_y,I_z] = size(Image_without_surf);
    Image_labels = zeros(I_x,I_y,I_z);
    for n = 1:I_z
        I_bin =  Image_without_surf(:,:,n);
        I2 = bwmorph(I_bin,'remove'); % deletes everything but surface voxels --> leaves pixels of new surface because outer surface was already removed!
        Image_labels(:,:,n) = I2(:,:);
        %         clear I_bin & I2 & n;
    end
    clear I_x & I_y & I_z;
    %% 2022-03-10 Update
    % We modified the code to calculate the surface at day 0 the same say they did for the other days, because D0 was
    % never done before.
    % change surface voxels day 0
    Result_volumes_inclsurf_d0_seg = (Result_volumes_inclsurf_d0 == 1) & (Bone_both == 1);
    % surface voxels:
    [I_x,I_y,I_z] = size(Result_volumes_inclsurf_d0_seg);
    Image_surf_labels_d0 = zeros(I_x,I_y,I_z);
    for n = 1:I_z %slicewise
        I_bin = Result_volumes_inclsurf_d0_seg(:,:,n);
        I2 = bwmorph(I_bin,'remove'); % deletes everything but surface voxels
        Image_surf_labels_d0(:,:,n) = I2(:,:);
        %         clear I_bin & I2 & n;
    end
    % cortex without surface for day 0
    
    Image_without_surf_logic_d0 = (Result_volumes_inclsurf_d0_seg == 1) & (Image_surf_labels_d0 == 0);
    Image_without_surf_d0 = uint8(Image_without_surf_logic_d0);
    [I_x,I_y,I_z] = size(Image_without_surf_d0);
    Image_labels_d0 = zeros(I_x,I_y,I_z);
    for n = 1:I_z
        I_bin =  Image_without_surf_d0(:,:,n);
        I2 = bwmorph(I_bin,'remove');
        Image_labels_d0(:,:,n) = I2(:,:);
        %         clear I_bin & I2 & n;
    end
    clear I_x & I_y & I_z;
    
    %% 2022-03-10 Update
    % Soft tissue interface
    Image_surf_const_labels_ST = Image_labels;
    % --------- NEW --------- %
    Image_surf_formed_labels_ST =  (Image_surf_const_labels_ST == 1) & (formed_bone == 1);
    Image_surf_formed_labels_ST = uint8(Image_surf_formed_labels_ST);
    
    %% 2022-03-10 Update from V1.5
    % Constent bone interface
    %  calculate Surface Area (SA) of resorbed:
    %  surface area of constant bone, grow 1, check if it´s newly formed
    [I_x,I_y,I_z] = size(Image_without_surf);
    Image_surf_const_labels_CB = zeros(I_x,I_y,I_z);
    Constant_bone_volume1 = (Image_without_surf == 1) & (constant_bone == 1);
    Constant_bone_volume = uint8(Constant_bone_volume1);
    for n = 1:I_z
        I_bin = Constant_bone_volume(:,:,n);
        disc_5 = strel('disk', 1,0);
        I2 = imdilate(I_bin, disc_5);
        I3 = bwmorph(I2, 'remove');
        Image_surf_const_labels_CB(:,:,n) = I3(:,:);
    end
    clear I_x & I_y & I_z;
    
    Image_surf_resorbed_labels_CB =  (Image_surf_const_labels_CB == 1) & (resorbed_bone == 1);
    Image_surf_resorbed_labels_CB = uint8(Image_surf_resorbed_labels_CB);
    
    % calculate SA of formed:
    Image_surf_formed_labels_CB =  (Image_surf_const_labels_CB == 1) & (formed_bone == 1);
    Image_surf_formed_labels_CB = uint8(Image_surf_formed_labels_CB);
    
    
    %% +++++++++++++++++get surface of formation/resorption sites++++++++++++++++
    % whole bone
    const_bone_surf_voxels = (Image_labels == 1) & (constant_bone == 1);
    formed_bone_surf_voxels_ST = Image_surf_formed_labels_ST;
    resorbed_bone_surf_voxels_ST = (Image_labels == 1) & (resorbed_bone == 1);
    %% 2022-03-10 Update
    formed_bone_surf_voxels_CB = Image_surf_formed_labels_CB;
    resorbed_bone_surf_voxels_CB = Image_surf_resorbed_labels_CB;
    %%
    switch which_eval
        case 'cort'
            % endosteal
            endostealArea = (Image_labels == 1) & (Endosteal_both == 1);
            endosteal_const_bone_surf_voxels = (endostealArea == 1) & (constant_bone == 1);
            endosteal_formed_bone_surf_voxels = (endostealArea == 1) & (formed_bone == 1);
            endosteal_resorbed_bone_surf_voxels = (endostealArea == 1) & (resorbed_bone == 1);
            
            % periosteal
            periostealArea = (Image_labels == 1) & (Periosteal_both == 1);
            periosteal_const_bone_surf_voxels = (periostealArea == 1) & (constant_bone == 1);
            periosteal_formed_bone_surf_voxels = (periostealArea == 1) & (formed_bone == 1);
            periosteal_resorbed_bone_surf_voxels = (periostealArea == 1) & (resorbed_bone == 1);
        case 'trab'
            % Do nothing
    end
    % get amount of surface voxels
    % whole bone
    const_bone_surf_voxels_nr = nnz(const_bone_surf_voxels);
    formed_bone_surf_voxels_nr_ST = nnz(formed_bone_surf_voxels_ST);
    resorbed_bone_surf_voxels_nr_ST = nnz(resorbed_bone_surf_voxels_ST);
    %% 2022-03-10 Update from V1.5
    formed_bone_surf_voxels_nr_CB = nnz(formed_bone_surf_voxels_CB);
    resorbed_bone_surf_voxels_nr_CB = nnz(resorbed_bone_surf_voxels_CB);
    %%
    total_bone_surf_voxels_nr = nnz(Image_labels_d0);
    switch which_eval
        case 'cort'
            % endosteal
            endosteal_const_bone_surf_voxels_nr = nnz (endosteal_const_bone_surf_voxels);
            endosteal_formed_bone_surf_voxels_nr = nnz (endosteal_formed_bone_surf_voxels);
            endosteal_resorbed_bone_surf_voxels_nr = nnz (endosteal_resorbed_bone_surf_voxels);
            % periosteal
            periosteal_const_bone_surf_voxels_nr = nnz (periosteal_const_bone_surf_voxels);
            periosteal_formed_bone_surf_voxels_nr = nnz (periosteal_formed_bone_surf_voxels);
            periosteal_resorbed_bone_surf_voxels_nr = nnz (periosteal_resorbed_bone_surf_voxels);
        case 'trab'
            % do nothing
    end
    % convert voxel into area: resolution voxelsize^2
    % whole bone
    % ---------- NEW ---------- %% 2019-12-13
    % we multiply by 1.5 now to compensate partial volume effect, first
    % layer was removed
    const_bone_surf_area = const_bone_surf_voxels_nr * (voxelsize^2) * 1.5;
    formed_bone_surf_area_ST = formed_bone_surf_voxels_nr_ST * (voxelsize^2) * 1.5;
    resorbed_bone_surf_area_ST = resorbed_bone_surf_voxels_nr_ST * (voxelsize^2) * 1.5;
    %% 2022-03-10 Update from V1.5
    formed_bone_surf_area_CB = formed_bone_surf_voxels_nr_CB * (voxelsize^2) * 1.5;
    resorbed_bone_surf_area_CB = resorbed_bone_surf_voxels_nr_CB * (voxelsize^2) * 1.5;
    %%
    total_bone_surf_area = total_bone_surf_voxels_nr * (voxelsize^2) * 1.5;
    switch which_eval
        case 'cort'
            % endosteal
            endosteal_const_bone_surf_area = endosteal_const_bone_surf_voxels_nr * (voxelsize^2) * 1.5;
            endosteal_formed_bone_surf_area = endosteal_formed_bone_surf_voxels_nr * (voxelsize^2) * 1.5;
            endosteal_resorbed_bone_surf_area = endosteal_resorbed_bone_surf_voxels_nr * (voxelsize^2) * 1.5;
            % periosteal
            periosteal_const_bone_surf_area = periosteal_const_bone_surf_voxels_nr * (voxelsize^2) * 1.5;
            periosteal_formed_bone_surf_area = periosteal_formed_bone_surf_voxels_nr * (voxelsize^2) * 1.5;
            periosteal_resorbed_bone_surf_area = periosteal_resorbed_bone_surf_voxels_nr * (voxelsize^2) * 1.5;
        case 'trab'
            % Do nothing
    end
    totalArea = total_bone_surf_area;
    %% 2022-03-10 Update
    % Added the _ST to formed and resorbed
    const_SATA = const_bone_surf_area/totalArea;
    formed_SATA_ST = formed_bone_surf_area_ST/totalArea;
    resorbed_SATA_ST = resorbed_bone_surf_area_ST/totalArea;
    %% 2022-03-10 Update 
    deltaArea_SATA_ST = formed_bone_surf_area_ST/totalArea + resorbed_bone_surf_area_ST/totalArea;

    %% 2022-03-10 Update from V1.5
    formed_SATA_CB = formed_bone_surf_area_CB/totalArea;
    resorbed_SATA_CB = resorbed_bone_surf_area_CB/totalArea;
    %% 2022-03-10 Update
    deltaArea_SATA_CB = formed_bone_surf_area_CB/totalArea + resorbed_bone_surf_area_CB/totalArea;

    %%
    clear totalArea;
    switch which_eval
        case 'cort'
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
        case 'trab'
            % Do nothing
    end
    %% calculate volumes
    % whole bone
    constant_bone_voxels_all_logic = (Image_without_surf==1) & (constant_bone==1);
    formed_bone_voxels_all_logic = (Image_without_surf==1) & (formed_bone==1);
    resorbed_bone_voxels_all_logic = (Image_without_surf==1) & (resorbed_bone==1);
    %% Creating new vis
    Result_volumes_inclsurf_labels = uint8(constant_bone_voxels_all_logic + ...
        formed_bone_voxels_all_logic*2 + ...
        resorbed_bone_voxels_all_logic*3);
    Result_volumes_inclsurf = uint8(constant_bone_voxels_all_logic + ...
        formed_bone_voxels_all_logic + ...
        resorbed_bone_voxels_all_logic);
    
    switch which_eval
        case 'cort'
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
        case 'trab'
    end
    % number of voxel in the total volume without the surface
    % whole bone
    formed_bone_voxels_all = nnz(formed_bone_voxels_all_logic);
    constant_bone_voxels_all = nnz(constant_bone_voxels_all_logic);
    resorbed_bone_voxels_all = nnz(resorbed_bone_voxels_all_logic);
    clear resorbed_bone_voxels_all_logic & constant_bone_voxels_all_logic & formed_bone_voxels_all_logic;
    
    switch which_eval
        case 'cort'
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
        case 'trab'
            % Do nothing
    end
    
    % number of voxel in the total volume with the surface -> half of
    % the surface voxel added due to partial volume effect -> the
    % original surface is not added!
    % whole bone
    constant_bone_voxels = constant_bone_voxels_all + const_bone_surf_voxels_nr*0.5;
    formed_bone_voxels = formed_bone_voxels_all + formed_bone_surf_voxels_nr_ST*0.5;
    resorbed_bone_voxels = resorbed_bone_voxels_all + resorbed_bone_surf_voxels_nr_ST*0.5;
    clear resorbed_bone_voxels_all & constant_bone_voxels_all & formed_bone_voxels_all;
    switch which_eval
        case 'cort'
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
        case 'trab'
            % Do nothing
    end
    % convert voxel into volume: resolution voxelsize^3
    % whole bone
    constant_bone_volume = constant_bone_voxels * voxelsize^3;
    formed_bone_volume = formed_bone_voxels * voxelsize^3;
    resorbed_bone_volume = resorbed_bone_voxels * voxelsize^3;
    clear resorbed_bone_voxels & formed_bone_voxels & constant_bone_voxels;
    switch which_eval
        case 'cort'
            % endosteal
            %     endosteal_constant_bone_volume = endosteal_constant_bone_voxels * voxelsize^3;
            endosteal_formed_bone_volume = endosteal_formed_bone_voxels * voxelsize^3;
            endosteal_resorbed_bone_volume = endosteal_resorbed_bone_voxels * voxelsize^3;
            clear endosteal_resorbed_bone_voxels & endosteal_formed_bone_voxels & endosteal_constant_bone_voxels
            % periosteal
            %     periosteal_constant_bone_volume = periosteal_constant_bone_voxels * voxelsize^3;
            periosteal_formed_bone_volume = periosteal_formed_bone_voxels * voxelsize^3;
            periosteal_resorbed_bone_volume = periosteal_resorbed_bone_voxels * voxelsize^3;
            clear periosteal_resorbed_bone_voxels & periosteal_formed_bone_voxels & periosteal_constant_bone_voxels
        case 'trab'
            % Do nothing
    end
    %VOLUME/total volume % total bone volume day 0
    % whole bone
    totalVolume     = constant_bone_volume + resorbed_bone_volume; % true volume at reference day!
    constant_BVTV   = constant_bone_volume / totalVolume;
    formed_BVTV  	= formed_bone_volume / totalVolume;
    resorbed_BVTV   = resorbed_bone_volume / totalVolume;
    %% 2022-03-10 Update
     deltaVolume = formed_BVTV + resorbed_BVTV;
     %%
    switch which_eval
        case 'cort'
            % endosteal
            endosteal_formed_BVTV = endosteal_formed_bone_volume / totalVolume;
            endosteal_resorbed_BVTV = endosteal_resorbed_bone_volume / totalVolume;
            % periosteal
            periosteal_formed_BVTV = periosteal_formed_bone_volume / totalVolume;
            periosteal_resorbed_BVTV = periosteal_resorbed_bone_volume / totalVolume;
            clear totalVolume;
        case 'trab'
            % Do nothing
    end
    % average height of volumes
    % whole bone
    formed_height = formed_bone_volume / formed_bone_surf_area_ST;
    resorb_height = resorbed_bone_volume / resorbed_bone_surf_area_ST;
    
    switch which_eval
        case 'cort'
            % endosteal
            endosteal_formed_height = endosteal_formed_bone_volume / endosteal_formed_bone_surf_area;
            endosteal_resorb_height = endosteal_resorbed_bone_volume / endosteal_resorbed_bone_surf_area;
            % periosteal
            periosteal_formed_height = periosteal_formed_bone_volume / periosteal_formed_bone_surf_area;
            periosteal_resorb_height = periosteal_resorbed_bone_volume / periosteal_resorbed_bone_surf_area;
        case 'trab'
            % Do nothing
    end
    
    % following the average thickness of formed and resorbed bone is calculated
    % mean thickness formed bone
    % whole bone
    I_formed_without_surf = double(Image_without_surf).* double(formed_bone);
    I_inv = -(I_formed_without_surf - 1);
    Dist = bwdist(I_inv,'euclidean'); % euclidian distance map of image
    Ult_E = bwulterode(Dist,'euclidean',6); % ultimate erosion of the distance map -> gets the local maxima -> point where most bone was formed
    E = double(Ult_E);
    dist_formed = 2 * double(Dist) - 1;
    loc_max_formed = E .* dist_formed * voxelsize; % Voxelsize
    % local maxima*value at local maxima = E*Dist; µm
    mean_formed =  mean (loc_max_formed(loc_max_formed~=0));
    std_formed = std (loc_max_formed(loc_max_formed~=0));
    max_formed = max(loc_max_formed(loc_max_formed~=0));
    clear I_formed_without_surf & I_inv & Dist & Ult_E & E;
    
    switch which_eval
        case 'cort'
            % endosteal
            I_formed_without_surf = double(Image_without_surf).* double(formed_bone) & (Endosteal_both == 1);
            I_inv = -(I_formed_without_surf - 1);
            Dist = bwdist(I_inv,'euclidean'); % euclidian distance map of image
            Ult_E = bwulterode(Dist,'euclidean',6); % ultimate erosion of the distance map -> gets the local maxima -> point where most bone was formed
            E = double(Ult_E);
            dist_formedEndosteal = 2 * double(Dist) - 1;
            loc_max_formedEndosteal = E .* dist_formedEndosteal * voxelsize;
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
            loc_max_formedPeriosteal = E .* dist_formedPeriosteal * voxelsize;
            % local maxima*value at local maxima = E*Dist; µm
            mean_formedPeriosteal = mean (loc_max_formedPeriosteal(loc_max_formedPeriosteal~=0));
            std_formedPeriosteal = std (loc_max_formedPeriosteal(loc_max_formedPeriosteal~=0));
            max_formedPeriosteal = max(loc_max_formedPeriosteal(loc_max_formedPeriosteal~=0));
            clear I_formed_without_surf & I_inv & Dist & Ult_E & E;
        case 'trab'
            % Do nothing
    end
    
    
    % mean thickness resorbed bone
    % whole bone
    I_resorbed_without_surf = double(Image_without_surf).* double(resorbed_bone);
    I_inv = -(I_resorbed_without_surf - 1);
    Dist = bwdist(I_inv,'euclidean'); % euclidian distance map of image
    Ult_E = bwulterode(Dist,'euclidean',6); % ultimate erosion of the distance map -> gets the local maxima -> point where most bone was formed
    E = double(Ult_E);
    dist_resorbed = 2 * double(Dist) - 1;
    loc_max_resorbed = E .* dist_resorbed * voxelsize;
    % local maxima*value at local maxima = E*Dist; µm
    mean_resorbed =  mean (loc_max_resorbed(loc_max_resorbed~=0));
    std_resorbed = std (loc_max_resorbed(loc_max_resorbed~=0));
    max_resorbed = max(loc_max_resorbed(loc_max_resorbed~=0));
    clear I_resorbed_without_surf & I_inv & Dist & Ult_E & E;
    
    switch which_eval
        case 'cort'
            % endosteal bone
            I_resorbed_without_surf = double(Image_without_surf).* double(resorbed_bone) & (Endosteal_both == 1);
            I_inv = -(I_resorbed_without_surf - 1);
            Dist = bwdist(I_inv,'euclidean'); % euclidian distance map of image
            Ult_E = bwulterode(Dist,'euclidean',6); % ultimate erosion of the distance map -> gets the local maxima -> point where most bone was formed
            E = double(Ult_E);
            dist_resorbedEndosteal = 2 * double(Dist) - 1;
            loc_max_resorbedEndosteal = E .* dist_resorbedEndosteal * voxelsize;
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
            loc_max_resorbedPeriosteal = E .* dist_resorbedPeriosteal * voxelsize;
            
            
            % local maxima*value at local maxima = E*Dist; µm
            mean_resorbedPeriosteal = mean (loc_max_resorbedPeriosteal(loc_max_resorbedPeriosteal~=0));
            std_resorbedPeriosteal = std (loc_max_resorbedPeriosteal(loc_max_resorbedPeriosteal~=0));
            max_resorbedPeriosteal = max(loc_max_resorbedPeriosteal(loc_max_resorbedPeriosteal~=0));
            clear I_resorbed_without_surf & I_inv & Dist & Ult_E & E;
        case 'trab'
            % Do nothing
    end
    
    % 3D - MAR
    day = cell2mat(days(m+1));
    % whole bone
    MAR_3D = mean_formed/day;
    MRR_3D = mean_resorbed/day;
    
    % Volume of bone formed/resorbed per time unit
    % whole bone
    BFR_3D = formed_bone_volume/day; % bone volume formed per day
    BRR_3D = resorbed_bone_volume/day;
    
    switch which_eval
        case 'cort'
            % endosteal
            endosteal_MAR_3D = mean_formedEndosteal / day;
            endosteal_MRR_3D = mean_resorbedEndosteal / day;
            % periosteal
            periosteal_MAR_3D = mean_formedPeriosteal / day;
            periosteal_MRR_3D = mean_resorbedPeriosteal / day;
            
            % Volume of bone formed/resorbed per time unit
            % endosteal
            endosteal_BFR_3D = endosteal_formed_bone_volume/day; % bone volume formed per day
            endosteal_BRR_3D = endosteal_resorbed_bone_volume/day;
            % periosteal
            periosteal_BFR_3D = periosteal_formed_bone_volume/day; % bone volume formed per day
            periosteal_BRR_3D = periosteal_resorbed_bone_volume/day;
        case 'trab'
            % Do nothing
    end
    
    
    
    %% Saving root
    switch which_eval
        case 'cort'
            path_bone_complete_VISres = strcat(char(res_folder(1)), '\cortical\');
        case 'trab'
            path_bone_complete_VISres = strcat(char(res_folder(1)), '\trabecular\');
            
    end
    
    if ~exist(path_bone_complete_VISres,'dir')
        mkdir(path_bone_complete_VISres);
    end
    % Write volumetric results
    filename = char(c_2(m));
    if is_short == 1
        %% 2022-03-12 Update
        % Add parameter for the two thresholds and the ST also
        switch which_eval
            case 'cort'
                Parameters = {'Threshold 0'; 'Threshold X'; 'R_1'; 'R_2'; 'R_3'; 'R_4'; 'Voxelsize'; 'Ec.MV/BV'; 'Ec.EV/BV'; 'Ec.MS/BS'; 'Ec.ES/BS';...
                    'Ec.MAR3D'; 'Ec.MRR3D'; 'Ec.BFR3D'; 'Ec.BRR3D'; 'Ec.MTh'; 'Ec.ED'; 'Ps.MV/BV'; 'Ps.EV/BV'; 'Ps.MS/BS'; 'Ps.ES/BS';...
                    'Ps.MAR3D'; 'Ps.MRR3D'; 'Ps.BFR3D'; 'Ps.BRR3D'; 'Ps.MTh'; 'Ps.ED'; 'T.MV/BV'; 'T.EV/BV'; 'T.MS_ST/BS' ; 'T.MS_CB/BS'; 'T.ES_ST/BS'; ...
                    'T.ES_CB/BS'; 'T.MAR3D'; 'T.MRR3D'; 'T.BFR3D'; 'T.BRR3D'; 'T.MTh'; 'T.ED'};
                
                Values = [level_0 ; level_x; R_1; R_2; R_3; R_4; voxelsize; endosteal_formed_BVTV; endosteal_resorbed_BVTV; endosteal_formed_SATA;...
                    endosteal_resorbed_SATA; endosteal_MAR_3D; endosteal_MRR_3D; endosteal_BFR_3D; endosteal_BRR_3D; ...
                    mean_formedEndosteal; mean_resorbedEndosteal; periosteal_formed_BVTV; periosteal_resorbed_BVTV; ...
                    periosteal_formed_SATA; periosteal_resorbed_SATA; periosteal_MAR_3D; periosteal_MRR_3D; periosteal_BFR_3D;...
                    periosteal_BRR_3D; mean_formedPeriosteal; mean_resorbedPeriosteal; formed_BVTV; resorbed_BVTV; formed_SATA_ST; formed_SATA_CB; ...
                    resorbed_SATA_ST; resorbed_SATA_CB; MAR_3D; MRR_3D; BFR_3D; BRR_3D; formed_height; resorb_height];
            case 'trab'
                Parameters = {'Threshold 0'; 'Threshold X'; 'R_1'; 'R_2'; 'R_3'; 'R_4'; 'Voxelsize'; 'T.MV/BV'; 'T.EV/BV'; 'T.MS_ST/BS' ; 'T.MS_CB/BS'; 'T.ES_ST/BS'; ...
                    'T.ES_CB/BS'; 'T.MAR3D'; 'T.MRR3D'; 'T.BFR3D'; 'T.BRR3D'; 'T.MTh'; 'T.ED'};
                
                Values = [level_0 ; level_x; R_1; R_2; R_3; R_4; voxelsize; formed_BVTV; resorbed_BVTV; formed_SATA_ST; formed_SATA_CB; ...
                    resorbed_SATA_ST; resorbed_SATA_CB; MAR_3D; MRR_3D; BFR_3D; BRR_3D; formed_height; resorb_height];
        end
        resultsfile = strcat(path_bone_complete_VISres,filename, '_',which_eval,'_VolResults_short.xlsx');
    else
        %% 2022-03-12 Update
        % Add parameter for the two thresholds and the ST also
        switch which_eval
            case 'cort'
                Parameters = {'Threshold 0'; 'Threshold X'; 'R_1'; 'R_2'; 'R_3'; 'R_4';'Voxelsize'; 'Ec.MV'; 'Ec.EV';...
                    'Ec.MV/BV'; 'Ec.EV/BV'; 'Ec.CS'; 'Ec.MS'; 'Ec.ES'; 'Ec.CS/BS'; 'Ec.MS/BS'; 'Ec.ES/BS'; ...
                    'Ec.MAR3D'; 'Ec.MRR3D'; 'Ec.BFR3D'; 'Ec.BRR3D'; 'Ec.MTh'; 'Ec.MTh.sd'; 'Ec.MTh.max'; ...
                    'Ec.ED'; 'Ec.ED.sd'; 'Ec.ED.max'; 'Ec.Md'; 'Ec.Ed';'Pt.MV'; 'Ps.EV' ; 'Ps.MV/BV';'Ps.EV/BV'; ...
                    'Ps.CS'; 'Ps.MS'; 'Ps.ES'; 'Ps.CS/BS'; 'Ps.MS/BS'; 'Ps.ES/BS'; 'Ps.MAR3D'; 'Ps.MRR3D'; 'Ps.BFR3D';...
                    'Ps.BRR3D'; 'Ps.MTh'; 'Ps.MTh.sd'; 'Ps.MTh.max'; 'Ps.ED'; 'Ps.ED.sd'; 'Ps.ED.max'; 'Ps.Md'; 'Ps.Ed';...
                    'T.CV'; 'T.MV'; 'T.EV'; 'T.CV/BV'; 'T.MV/BV'; 'T.EV/BV'; 'T.CS'; 'T.MS'; 'T.ES'; 'T.CS/BS'; 'T.MS_ST/BS'; 'T.MS_CB/BS'; 'T.ES_ST/BS';...
                    'T.ES_CB/BS'; 'T.MAR3D'; 'T.MRR3D'; 'T.BFR3D'; 'T.BRR3D'; 'T.MTh'; 'T.MTh.sd'; 'T.MTh.max'; 'T.ED'; 'T.ED.sd';...
                    'T.ED.max'; 'T.Md'; 'T.Ed';'DeltaBV'; 'DeltaBS_ST';'DeltaBS_CB'};
                Values = [level_0; level_x; R_1; R_2; R_3; R_4; voxelsize; endosteal_formed_bone_volume;...
                    endosteal_resorbed_bone_volume; endosteal_formed_BVTV; endosteal_resorbed_BVTV; ...
                    endosteal_const_bone_surf_area; endosteal_formed_bone_surf_area; endosteal_resorbed_bone_surf_area; ...
                    endosteal_const_SATA; endosteal_formed_SATA; endosteal_resorbed_SATA; endosteal_MAR_3D; endosteal_MRR_3D;...
                    endosteal_BFR_3D; endosteal_BRR_3D; mean_formedEndosteal; std_formedEndosteal; max_formedEndosteal; ...
                    mean_resorbedEndosteal; std_resorbedEndosteal; max_resorbedEndosteal; endosteal_formed_height; endosteal_resorb_height;...
                    periosteal_formed_bone_volume; periosteal_resorbed_bone_volume; periosteal_formed_BVTV; periosteal_resorbed_BVTV; ...
                    periosteal_const_bone_surf_area; periosteal_formed_bone_surf_area; periosteal_resorbed_bone_surf_area; ...
                    periosteal_const_SATA; periosteal_formed_SATA; periosteal_resorbed_SATA; periosteal_MAR_3D; periosteal_MRR_3D; ...
                    periosteal_BFR_3D; periosteal_BRR_3D; periosteal_formed_height; periosteal_resorb_height; mean_formedPeriosteal;...
                    std_formedPeriosteal; max_formedPeriosteal; mean_resorbedPeriosteal; std_resorbedPeriosteal; max_resorbedPeriosteal;...
                    constant_bone_volume; formed_bone_volume; resorbed_bone_volume; constant_BVTV; formed_BVTV; resorbed_BVTV; ...
                    const_bone_surf_area; formed_bone_surf_area_ST; resorbed_bone_surf_area_ST; const_SATA; formed_SATA_ST; formed_SATA_CB; ...
                    resorbed_SATA_ST; resorbed_SATA_CB; MAR_3D; MRR_3D; BFR_3D; BRR_3D; mean_formed; std_formed; max_formed; mean_resorbed; ...
                    std_resorbed; max_resorbed; formed_height; resorb_height; deltaVolume; deltaArea_SATA_ST; deltaArea_SATA_CB];
            case 'trab'
                Parameters = {'Threshold 0'; 'Threshold X'; 'R_1'; 'R_2'; 'R_3'; 'R_4'; 'Voxelsize';...
                    'T.CV'; 'T.MV'; 'T.EV'; 'T.CV/BV'; 'T.MV/BV'; 'T.EV/BV'; 'T.CS'; 'T.MS'; 'T.ES'; 'T.CS/BS'; 'T.MS_ST/BS' ; 'T.MS_CB/BS'; 'T.ES_ST/BS'; ...
                    'T.ES_CB/BS'; 'T.MAR3D'; 'T.MRR3D'; 'T.BFR3D'; 'T.BRR3D'; 'T.MTh'; 'T.MTh.sd'; 'T.MTh.max'; 'T.ED'; 'T.ED.sd';...
                    'T.ED.max'; 'T.Md'; 'T.Ed'; 'DeltaBV'; 'DeltaBS_ST';'DeltaBS_CB'};
                Values = [level_0; level_x; R_1; R_2; R_3; R_4; voxelsize; ...
                    constant_bone_volume; formed_bone_volume; resorbed_bone_volume; constant_BVTV; formed_BVTV; resorbed_BVTV; ...
                    const_bone_surf_area; formed_bone_surf_area_ST; resorbed_bone_surf_area_ST; const_SATA; formed_SATA_ST; formed_SATA_CB; ...
                    resorbed_SATA_ST; resorbed_SATA_CB; MAR_3D; MRR_3D; BFR_3D; BRR_3D; mean_formed; std_formed; max_formed; mean_resorbed; ...
                    std_resorbed; max_resorbed; formed_height; resorb_height; deltaVolume; deltaArea_SATA_ST; deltaArea_SATA_CB];
        end
        
        resultsfile = strcat(path_bone_complete_VISres, filename, '_',which_eval,'_VolResults_long.xlsx');
    end
    T = table(Parameters, Values);
    writetable(T, resultsfile);
    
    %% Visuaslization step (create the outputs)
    switch which_eval
        case 'cort'
        case 'trab'
            % create fake endo and peri files because they dont
            % exist in trabecular but the functions need an
            % input, that will be empty cell.
            Endosteal_both = [];
            Periosteal_both  = [];
    end
    
    if is_short == 1 % if the report is short
        if is_mat == 1 % If you ask for saving in .mat file
            %% THIS HAS NOT BEEN UPDATED, SO MAY NOT BE COMPLETE.
            warning('Saving mat file has not been updated, so .mat file may be incomplete');
            %% Old
            file_bone_complete_VISres = strcat(filename, '_',which_eval,'_complete_ImageResults_short');
            savefile2 = ([path_bone_complete_VISres file_bone_complete_VISres]);
            save(savefile2, 'Result_volumes_inclsurf_labels' , 'Result_volumes_inclsurf');
        else % is tif
            if strcmp(VisualizationSwitch, 'Short') % If the outputs needed are the short, irrespective of the full/short report
                save_visualization_short(day, path_bone_complete_VISres, Result_volumes_inclsurf_labels, Result_volumes_inclsurf);
            elseif strcmp(VisualizationSwitch , 'Full')
                
                save_visualization_full(which_eval, day, path_bone_complete_VISres, Bone_both, Image_t0, Image_tX, ...
                    Result_volumes_inclsurf_labels, Result_volumes_inclsurf, constant_bone, formed_bone, resorbed_bone, ...
                    Endosteal_both, Periosteal_both, dist_formed, loc_max_formed, dist_resorbed, loc_max_resorbed)
            end
        end
    else
        if is_mat == 1
            %% THIS HAS NOT BEEN UPDATED, SO MAY NOT BE COMPLETE.
            warning('Saving mat file has not been updated, so .mat file may be incomplete');
            %% Old
            file_bone_complete_VISres = strcat(filename, '_',which_eval,'_complete_ImageResults_long');
            savefile2 = ([path_bone_complete_VISres file_bone_complete_VISres]);
            save(savefile2, 'Bone_both', 'Image_t0', 'Image_tX', 'Result_volumes_inclsurf_labels' , 'Result_volumes_inclsurf');
            save(savefile2,  'Image_without_surf','-append' )
            save(savefile2,'constant_bone','formed_bone','resorbed_bone','-append' );
            save(savefile2, 'Endosteal_both' ,'Periosteal_both','-append')
            save(savefile2,'dist_formed','loc_max_formed','dist_resorbed','loc_max_resorbed','-append' );
        else % is tif
            if strcmp(VisualizationSwitch, 'Short') % If the outputs needed are the short, irrespective of the full/short report
                save_visualization_short(day, path_bone_complete_VISres, Result_volumes_inclsurf_labels, Result_volumes_inclsurf);
            elseif strcmp(VisualizationSwitch , 'Full')
                save_visualization_full(which_eval, day, path_bone_complete_VISres, Bone_both, Image_t0, Image_tX, ...
                    Result_volumes_inclsurf_labels, Result_volumes_inclsurf, constant_bone, formed_bone, resorbed_bone, ...
                    Endosteal_both, Periosteal_both, dist_formed, loc_max_formed, dist_resorbed, loc_max_resorbed)
            end
        end
    end
end
%% Ending
close(progress);
close(uf);
toc
end
