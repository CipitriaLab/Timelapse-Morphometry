function folder = imwrite3D_stacked(I,varargin)
%% imwrite3D_stacked(I,varargin)
% I is the 3D image
% First arg: is the folder
% Second arg: is the file
% third is the voxelsize
%% 2021-05-07 Update
% Now the fourth argument can be voxelsize to store it for amira format, so
% I also updated the way to save
%%
if isempty(varargin)
    [~, folder] = get_old_folder_position();
elseif numel(varargin) == 1
    folder = varargin{1,1};
    if isempty(folder)
        [~, folder] = get_old_folder_position();
    end
    base_name = ask_input_filename();
elseif numel(varargin) == 2
    folder = varargin{1,1};
    if isempty(folder)
        [~, folder] = get_old_folder_position();
    end
    base_name = varargin{2};
elseif numel(varargin) == 3
    folder = varargin{1,1};
    if isempty(folder)
        [~, folder] = get_old_folder_position();
    end
    base_name = varargin{2};
    if isempty(base_name)
        base_name = ask_input_filename();
    end
    % then the voxelsize was input
    voxelsize = varargin{3};
end

%folder = path2clean(folder);
% filename = strcat(folder,base_name,'_','stacked'); % only in general
filename = strcat(folder,base_name); % not stacked for 3D timelapse writting
% 2021-07-13 Update
% if working on network drives, it is slow to access all the time to
% rewrite, so write locally and then copy. After a quick test it is so so
% so much faster it is incredible!
if isgnu() == 0
    fast_tmp = path2clean('C:/tmp/');
else
    fast_tmp = path2clean('/scratch/tmp/');
end
if exist(fast_tmp,'dir') == 0
    mkdir(fast_tmp);
end
filename_real = filename;
filename = strcat(fast_tmp, base_name);

msg = char(strcat('Creating 3D stack:'," ",base_name,'.tif'));
disp(msg);
comp_method = 'lzw'; % lossless method, so no problem for data
%% Depending on the input
if isa(I,'cell')
    % Watch out, the stacking order matter
    if size(I,1) == 1
        nb_frame = size(I,2);
        I = I';
    else
        nb_frame = size(I,1);
    end
    I_to_write = I;
else
    if numel(size(I)) == 4 % if it is RGB STACKS
        nb_frame = size(I,4);
        I_to_write = cell(nb_frame,1);
        for k = 1 : nb_frame
            I_to_write{k,1} = I(:,:,:,k);
        end
    else
        nb_frame = size(I,3);
        I_to_write = cell(nb_frame,1);
        for k = 1 : nb_frame % GRAYSCALE STACK
            I_to_write{k,1} = I(:,:,k);
        end
    end
end
%% Do the actual writting
k = 1;
dim = [size(I_to_write{k,1},1), size(I_to_write{k,1},2), size(I_to_write{k,1},3)]; % does not matter for stack RGB
for k = 1 : nb_frame
    if k == 1
        if exist('voxelsize','var')
            vox_x = voxelsize(1);
            vox_y = voxelsize(2);
            vox_z = voxelsize(3);
            xlim = round(vox_x * (dim(1) - 1),4);
            ylim = round(vox_y * (dim(2) - 1),4);
            zlim = round(vox_z * (nb_frame - 1),4);
%             disp(xlim); disp(ylim); disp(zlim);
            %% Watch out, from amira the x and y are switched! 
            descript = char(strcat('"BoundingBox'," ",'0'," ",num2str(ylim)," ",'0'," ",num2str(xlim)," ",'0'," ",num2str(zlim),'"'));
            disp(descript)
            imwrite(I_to_write{k,1}, strcat(filename,'.tif'),'Compression',comp_method,'Description',descript);
        else
            imwrite(I_to_write{k,1}, strcat(filename,'.tif'),'Compression',comp_method);
        end
    else
        imwrite(I_to_write{k,1},  strcat(filename,'.tif'),'WriteMode','append', 'Compression',comp_method);
    end
    display_progress(nb_frame,4,k,'Saving 3D stack: ');
end
movefile(strcat(filename,'.tif'), strcat(filename_real,'.tif'));
% explorer_smarter(folder);
end
