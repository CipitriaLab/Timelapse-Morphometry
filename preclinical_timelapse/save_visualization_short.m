function save_visualization_short(day, str_path, Result_volumes_inclsurf_labels, Result_volumes_inclsurf)
% ======================================================================= %
%% 2022-03-10
% save_visualization_short(day, str_path, Result_volumes_inclsurf_labels, Result_volumes_inclsurf)
% Will create the short visualization outputs (2 outputs).
% ======================================================================= %
Day_cell = num2str(day); % get the day number for the filename

str_file = strcat('Result_volumes_inclsurf_labels_', Day_cell); % create the base filename
imwrite3D_stacked(Result_volumes_inclsurf_labels, str_path, str_file);  % save the stack

str_file_RVIS = strcat('Result_volumes_inclsurf_', Day_cell);
imwrite3D_stacked(Result_volumes_inclsurf, str_path, str_file_RVIS);

end