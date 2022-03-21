function[] = static(inputfiles, paths, thresh, morph, bone_main_axis, voxelsize, res_folder)
%% 2022-03-10 Update
% This comes V1.6 "seg_only.m"
% ======================================================================= %
% INPUT
% inputfiles - nx1 cell array, ordered according to chronology of n scans
% paths - nx1 cell array, ordered as files
% thresh - cortical threshold, double
% morph - 4x1 double array with morphometric radii

% open up a progress bar to check whether the programm is running
uf = uifigure;
uf.Name = 'Running Evaluation';
uf.InnerPosition = [500 500 400 75];
progress = uiprogressdlg(uf,'Message','Please wait...','Indeterminate','on');

% Get files
reffile = inputfiles{1};
path = paths{1};
%% Change Hub and Sarah
level_x = thresh(2,1);

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


%% 
if b_1 == 1
    IMG_dayX = importdata([path_2 char(files(1))]);
    Image_dayX(:,:,:) = IMG_dayX(1,:,:,:);
else % input is tiff, slice-wise opening
    IMG_dayX = Tiff([path_2 char(files(1))], 'r');
    infox = imfinfo([path_2 char(files(1))]);
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
 Image_dayX = switching_direction(Image_dayX, bone_main_axis);

%% NO SMOOTH
% ------ NO SMOOTH -------% 
Image_X = Image_dayX;
% ------ NO SMOOTH -------% 
Image_treshX = Image_X >= level_x;
clear IMG_dayX & Image_dayX ; 

% Convert file X to 8-bit
Image_tX = uint8(Image_treshX);
% clear Image_treshX;

Result_volumes_inclsurf = uint8(Image_tX);
% create result arrays and fill them with zeros
[I_x,I_y,I_z] = size(Result_volumes_inclsurf);
Result_trab =   zeros(I_x,I_y,I_z);
Result_cort =   zeros(I_x,I_y,I_z);
Endosteal_both = zeros(I_x,I_y,I_z);
Periosteal_both = zeros(I_x,I_y,I_z);

for n = 1:I_z
    I_bin = Result_volumes_inclsurf(:,:,n); % this is the binary image of the bone (slicewise)
    %     dilation followed by closing of the image, smooth it a bit more and close little holes or
    %     lesions propagating through the cortex
    I2 = imdilate(I_bin,disc_1);
    I3 = imclose(I2,disc_2);
    %     removing holes in the images, -> filling the medullary canal
    I4 = imfill(I3,8,'holes');
    I5 = imerode(I4,disc_3);
    %     compare the filled binary image with the original one, this seperates
    %     the medullary canal from the cortex
    I6 = (I_bin == 0) & (I5 == 1);
    %     & creates and logical array, change it back to uint8
    I7 = uint8(I6);
    %     another closing operation is needed to fill the gaps in the
    %     cancellous bone mask
    I8 = imclose(I7,disc_4);
    %     erode again to remove any cortical bone from the cancellous bone mask
    I9 = imerode(I8,disc_1);
    %     comparison between the original binary image and the cancellous bone
    %     mask to seperate trabecular and cortical bone
    I10 = (I_bin == 1) & (I9 == 1);
    trab_bone = uint8(I10);
    I11 = (I_bin == 1) & (I9 == 0);
    cort_bone = uint8(I11);
    % Separating endosteum and periosteum
    % dilate first the cortical bone mask
    cort_filled = imdilate(cort_bone, disc_4);
    % filling the medullary canal
    I12 = imfill(cort_filled,'holes');
    % erode again with same disc
    cort_filled2 = imerode(I12, disc_4);
    % use shrink operation to shrink mask
    I13 = bwmorph(cort_filled2,'shrink',3);
    % Endosteum
    I14 = (cort_bone == 1) & (I13 == 1);
    % Periosteum
    I15 = (cort_bone == 1) & (I13 == 0);
    I16 = uint8(I14);
    I17 = uint8(I15);
    Result_TV_trab(:,:,n) = I9(:,:);
    Result_trab(:,:,n) = trab_bone(:,:);
    Result_TV_cort(:,:,n) = I4(:,:);
    Result_cort(:,:,n) = cort_bone(:,:);
    Endosteal_both(:,:,n) = I16(:,:);
    Periosteal_both(:,:,n) = I17(:,:);
end

C = zeros(I_x,I_y,I_z);
D = zeros(I_x,I_y,I_z);
for n = 1:I_z
    C(:,:,n) = Result_cort(:,:,n) + 2*Result_trab(:,:,n);
    D(:,:,n) = Endosteal_both(:,:,n) + 2*Periosteal_both(:,:,n);
end

%% surface calcs
% here follows the surface calculations

[I_x,I_y,I_z] = size(Result_cort);
    Image_surf_labels = zeros(I_x,I_y,I_z);
    parfor n = 1:I_z %slicewise
        I_bin = Result_cort(:,:,n);
        I2 = bwmorph(I_bin,'remove'); % deletes everything but surface voxels
        Image_surf_labels(:,:,n) = I2(:,:);
    end
    clear I_x & I_y & I_z;
    cort_BS = nnz(Image_surf_labels) *(voxelsize^2);

%% CALCULATIONS
% here follows calculations for BV/TV and BS/BV for both cortical bone and
% trabecular bone
trab_TV = nnz(Result_TV_trab) * (voxelsize^3);
trab_BV = nnz(Result_trab) * (voxelsize^3);
trab_BVTV = trab_BV/trab_TV;
%cortical bone
cort_TV = nnz(Result_TV_cort) * (voxelsize^3);
cort_BV = nnz(Result_cort) * (voxelsize^3);
cort_BVTV = cort_BV/cort_TV;
cort_BSBV = cort_BS/cort_BV;

%% write data
path_VISres = strcat(char(res_folder(1)), '\static\');
if ~ exist(path_VISres,'dir')
    mkdir(path_VISres);
end
filename = char(c_2(1));
Parameters = {'trab_TV';'trab_BV';'trab_BVTV';'cort_TV';'cort_BV';'cort_BVTV';'cort_BSBV'};
Values = [trab_TV;trab_BV;trab_BVTV;cort_TV;cort_BV;cort_BVTV;cort_BSBV];
resultsfile = strcat(path_VISres, filename, '_static.xlsx');

T = table(Parameters, Values);
writetable(T, resultsfile);

% close the progress bar
close(progress);
close(uf);
end