function folder_correct = path2clean(folder)
% This function is to make sure the folder path are always compatible with
% GNU/Linux and also that they end with the separateor
% 2020-05-23 
% Now support also if you specify a full filename and not just a path.

if contains(folder,'.')
    % New approach, just use explorer and try to open the folder to see
    % what is here
    if isfolder(folder)
        %     idx_dot = find(folder == '.');
        %     if idx_dot < numel(folder) - 10
        %         % then it is not the extension
    else
        if isfile(folder)
            % In this case folder is actually a full filename.
            [pathy,file,ext] = fileparts(folder);
            pathy = path2clean(pathy);
            folder_correct = strcat(pathy,file,ext);
%             disp('You specify a filename in path2clean instead of a folder');
            return;
        else
            msg = char(strcat(folder," ",'doest not exist anymore'));
            disp(msg);
        end
    end
end

folder(folder == '\') = '/';
if folder(end) ~= '/'
    folder_correct = strcat(folder,'/');
else
    folder_correct = folder;
end


end