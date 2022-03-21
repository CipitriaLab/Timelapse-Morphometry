% RUNS IN MATLAB VERSION R2016b 9.1.0.441655
% it creates an interface where you can select your structure element
% parameters
% calls an implay function at the end to look slicewise through the
% segmentation

function[] = seg_only(inputfiles, paths, thresh, morph, bone_main_axis)
% INPUT
% inputfiles - nx1 cell array, ordered according to chronology of n scans
% paths - nx1 cell array, ordered as files
% thresh - cortical threshold, double
% morph - 4x1 double array with morphometric radii

% open up a progress bar to check whether the programm is running
uf = uifigure;
name = uf.Name;
uf.Name = 'Running Evaluation';
uf.InnerPosition = [500 500 400 75];
progress = uiprogressdlg(uf,'Message','Please wait...','Indeterminate','on');

% Get files
reffile = inputfiles{1};
path = paths{1};
%% Change Hub and Sarah
% level = thresh; % old 
%level_0 = thresh(1,1);
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

% Open ref file
%% TO UNCOMMENT TO HAVE THE REF
% % if b_1 == 1
% %     IMG_day0 = importdata([path reffile]);
% %     Image_day0(:,:,:) = IMG_day0(1,:,:,:);
% % else % input is tiff, slice-wise opening
% %     IMG_day0 = Tiff([path reffile], 'r');
% %     info = imfinfo([path reffile]);
% %     s = size(info);
% %     s = s(1);
% %     r1 = read(IMG_day0);
% %     S = size(r1);
% %     m = S(1);
% %     n = S(2);
% %     Image_day0 = zeros(m, n, s);
% %     for i = 1:s-1
% %         Image_day0(:,:,i) = read(IMG_day0);
% %         nextDirectory(IMG_day0);
% %     end
% %     Image_day0 = int16(Image_day0);
% % end
% % % switching x,y,z directions to y,z,x -> different position for
% % % segmentation
% % Image_day0 = switching_direction(Image_day0, bone_main_axis);
% % 
% % % Filtering ref file
% % Image_0 = smooth3(Image_day0,'gaussian',[3 3 3],0.65);
% % % Image_0 = Image_day0;
% % Image_tresh0 = Image_0 >= level_0;
% % clear IMG_day0 & Image_day0 & Image_0;
% % 
% % % Convert ref file into 8bit integer
% % Image_t0 = uint8(Image_tresh0);
% % clear Image_tresh0;

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
  %  Image_dayX = int16(Image_dayX);
% Image_dayX = Image_dayX;
end
 Image_dayX = switching_direction(Image_dayX, bone_main_axis);
% min_to_set  = min(min(min(Image_dayX)));
% Image_dayX = uint16(Image_dayX + abs(min_to_set));
  
%  R_1 = Image_dayX(:,n,:)
% Filtering/thresholding file X
%%  SMOOTH
% ------ SMOOTH -------% 
% Image_X = smooth3(Image_dayX,'gaussian',[3 3 3],0.65);
% ------ SMOOTH -------% 
%% NO SMOOTH
% ------ NO SMOOTH -------% 
Image_X = Image_dayX;
% ------ NO SMOOTH -------% 
Image_treshX = Image_X >= level_x;
clear IMG_dayX & Image_dayX ; 

% Convert file X to 8-bit
Image_tX = uint8(Image_treshX);
% clear Image_treshX;

%% Goodbye 
% resorbed bone = volumes only in day 0
% resorbed_bone_X = (Image_tX == 0) & (Image_t0 == 1);
% resorbed_bone = uint8(resorbed_bone_X);
% 
% % formed bone = volumes only in day X
% formed_bone_X = (Image_tX == 1) & (Image_t0 == 0);
% formed_bone = uint8(formed_bone_X);
% 
% % constant bone = volumes in both datasets
% constant_bone_X = Image_t0 - resorbed_bone;
% constant_bone = uint8(constant_bone_X);

% Creating labels
% Result_volumes_inclsurf = constant_bone + formed_bone + resorbed_bone;
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
    Result_trab(:,:,n) = trab_bone(:,:);
    Result_cort(:,:,n) = cort_bone(:,:);
    Endosteal_both(:,:,n) = I16(:,:);
    Periosteal_both(:,:,n) = I17(:,:);
    %     clear I_bin & I2 & I3 & I4 & I5 & I6 & I7 & I8 & I9 & I10 & I11 & I12 & I13 & I14 & I15 & I16 & I17 & trab_bone & cort_bone
end

Result_trab_bin = logical(Result_trab);
Result_trab_8bits = uint8(Result_trab);
% to_write = Image_X;
% to_write(Result_trab_bin == 0) = 0;
% to_write_str = strcat(c_1,'_seg');

%% Mask , not neeeded anymore
% I_norm = double(to_write)./2^16;
% I_8bits = uint8(I_norm.*255);
% folder_str = imwrite3D_seq(I_8bits,[],to_write_str);

%% Saving the trab binary
folder_root = path2clean(paths{1});
to_write_bin_str = strcat(c_1, '_seg-binary-trab');
folder_str = strcat(folder_root,'seg_binary_trab/');
if exist(folder_str, 'dir') == 0
    mkdir(folder_str);
end
imwrite3D_seq(Result_trab_8bits,folder_str,to_write_bin_str);

%% Saving the cort binary
Result_cort_bin = logical(Result_cort);
Result_cort_8bits = uint8(Result_cort);
to_write_bin_cort_str = strcat(c_1, '_seg-binary-cort');
folder_str = strcat(folder_root,'seg_binary_cort/');
if exist(folder_str, 'dir') == 0
    mkdir(folder_str);
end
imwrite3D_seq(Result_cort_8bits,folder_str,to_write_bin_cort_str);

% % the following is just to check if the image segmentation went well!

C = zeros(I_x,I_y,I_z);
D = zeros(I_x,I_y,I_z);
for n = 1:I_z
    C(:,:,n) = Result_cort(:,:,n) + 2*Result_trab(:,:,n);
    D(:,:,n) = Endosteal_both(:,:,n) + 2*Periosteal_both(:,:,n);
end



% close the progress bar
close(progress);
close(uf);

% % define the limits for the colormap in implay
% limits=[0 2];
% ImplayWithMap(C, 50, limits);
% ImplayWithMap(D, 50, limits);
% 
% 
% 
%     function [handle] = ImplayWithMap(frames, fps, limits)
%         %ImplayWithMap Calls the implay function and adjust the color map
%         % Call it with 3 parameters:
%         % ImplayWithMap(frames, fps, limits)
%         % frames - 4D arrray of images
%         % fps - frame rate
%         % limits - an array of 2 elements, specifying the lower / upper
%         % of the liearly mapped colormap
%         % Returns a nadle to the player
%         %
%         % example:
%         % h = ImplayWithMap(MyFrames, 30, [10 50])
%         
%         
%         handle = implay(frames, fps);
%         handle.Visual.ColorMap.UserRangeMin = limits(1);
%         handle.Visual.ColorMap.UserRangeMax = limits(2);
%         handle.Visual.ColorMap.UserRange = 1;
%     end

end