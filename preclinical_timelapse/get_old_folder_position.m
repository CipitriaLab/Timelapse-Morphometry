function [filename,path_folder] = get_old_folder_position(varargin)

% This function checks if the file "defaut_path.mat" exist in the current
% directory. If it is it and the day stored is the current day, it will load the previous position in the tree
% directory of windows explorer.

% Argument : 1 ) extension is a string, if nothing is put, the uigetdir will be
% launched, meaning selection of a folder.
%           example of extension : 'tif', 'png', 'jpeg'
%% Update 2019-01-27 Fixing the overwritting and last extension used
if nargin == 0
    extension = [];
else
    extension = varargin;
end

c = clock ; % Get the current time
date_today = c(1:3);
ex_def_path = exist(strcat(pwd,'\default_path.mat'),'file');
%%%%%%%%%%%%%%%%%%%%
%% update 2019-01-27 I comment this line which was useless
% date_old = date_today; %
%%%%%%%%%%%%%%%%%%
%% If extension = [] ;
if isempty(extension)
    filename = [];
    switch ex_def_path
        case 0
            path_folder = uigetdir;
            if isnumeric(path_folder)
                error('The user stop the selection');
            end
            % If the user stopped the selection, an error message will pop
            % and nothing will be saved!
            path_folder_saved = path_gnu_win(path_folder);
            date_saved = date_today;
            extension_saved = extension;
            save('default_path.mat','path_folder_saved','date_saved','extension_saved');
        case 2
            load('default_path.mat','path_folder_saved','date_saved','extension_saved');
            test_data = date_saved == date_today;
            switch sum(test_data)
                case 3 % meaning the year, month and date are the same
                    path_folder = uigetdir(path_folder_saved,'Select a folder');
                    if isnumeric(path_folder) % meaning if the user did not cancel the selection
                        error('The user stopped the selection')
                    end % if the user pressed cancel, path_folder will be egal to 0
                    
                otherwise % The file exist but the path is from another day, start from root then
                    path_folder = uigetdir;
                    if isnumeric(path_folder)
                        error('The user stop the selection');
                    end
            end
            path_folder_saved = path_gnu_win(path_folder);
            date_saved = date_today;
            extension_saved = extension;
            save('default_path.mat','path_folder_saved','date_saved','extension_saved');
        otherwise
            % super unlikly to enter this problem
            error('The path specified is wrong');
    end
    %% Else the extension is specified
else
    extension_full = cell(nargin,1);
    for ee = 1 : nargin
        extension_full(ee,1) = strcat('*.',extension(ee));
    end
    switch ex_def_path
        case 0
            % Take care, the uigetfile add the last slash !
            [filename,path_folder] = uigetfile(extension_full,'File Selector'); % open a window and select your image
            if isnumeric(path_folder)
                error('The user stop the selection');
            end
            if isgnu()
                filename_full_saved  = strcat(path_folder,filename);
                filename_full_saved(filename_full_saved =='\') = '/';
            else
                filename_full_saved  = strcat(path_folder,filename);
            end
            date_saved = date_today;
            idx_extension = find(filename_full_saved =='.');
            if isempty(idx_extension)
                error('Did not find the extension');
            end
            extension_saved = filename_full_saved(idx_extension(end)+1:end);
            % by tacking idx_extension(end), one take care that even if
            % there are several points in the filename, the last will be
            % used
            path_folder_saved = path_gnu_win(path_folder); % to add the last slack/blackslash to the folder
            save('default_path.mat','path_folder_saved','date_saved','extension_saved','filename_full_saved');
        case 2
            tmp = load('default_path.mat');
            if isfield(tmp,'filename_full_saved')
                clear tmp
                load('default_path.mat','path_folder_saved','date_saved','extension_saved','filename_full_saved');
            else
                load('default_path.mat','path_folder_saved','date_saved','extension_saved');
            end
            test_date = date_saved == date_today ;
            switch sum(test_date)
                case 3
                    idx_match = find(strcmp(extension_saved,extension) ==1);
                    if ~isempty(idx_match) &&  numel(extension_saved) == numel(extension{1,idx_match})
                        % for tif and tiff for example, it really has to be the same extension to go here
                        [filename,path_folder] = uigetfile(filename_full_saved,'File Selector');
                        if isnumeric(path_folder)
                            error('The user stop the selection');
                        end
                    else
                        [filename,path_folder] = uigetfile(strcat(path_folder_saved,extension_full{1,1}),'File Selector');
                        if isnumeric(path_folder)
                            error('The user stop the selection');
                        end
                    end
                otherwise
                    [filename,path_folder] = uigetfile(strcat(path_folder_saved,extension_full{1,1}),'File Selector');
                    if isnumeric(path_folder)
                        error('The user stop the selection');
                    end
            end
             if isgnu()
                filename_full_saved  = strcat(path_folder,filename);
                filename_full_saved(filename_full_saved =='\') = '/';
            else
                filename_full_saved  = strcat(path_folder,filename);
            end
            date_saved = date_today;
            idx_extension = find(filename_full_saved =='.');
            if isempty(idx_extension)
                error('Did not find the extension');
            end
            extension_saved = filename_full_saved(idx_extension(end)+1:end);
            % by tacking idx_extension(end), one take care that even if
            % there are several points in the filename, the last will be
            % used
            path_folder_saved = path_gnu_win(path_folder); % to add the last slack/blackslash to the folder
            save('default_path.mat','path_folder_saved','date_saved','extension_saved','filename_full_saved');
        otherwise
            error('The path specified is wrong');
            
    end
    
end

end


function path_correct = path_gnu_win(path_to_correct)

if isgnu()
    path_correct = strcat(path_to_correct,'/');
else
    path_correct = strcat(path_to_correct,'\');
end
end

function out = isgnu()

path_test = pwd;
if strcmp(path_test(1:5),'/usr/')
    out = true;
else 
    out = false;
end


end

