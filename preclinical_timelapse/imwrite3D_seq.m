function folder = imwrite3D_seq(I,varargin)

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
end

folder = path2clean(folder);
% base_name = ask_input_filename();
for k = 1 : size(I,3)
   filename = strcat(folder,base_name,'_',get_digit(k));
   imwrite(I(:,:,k), strcat(filename,'.tif'));
   display_progress(size(I,3),4,k,'Saving 3D sequence: ');
end

end