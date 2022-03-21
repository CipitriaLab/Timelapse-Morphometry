function Image = switching_direction(Image, bone_main_axis)
% ======================================================================= %
%% 2022-03-10 Update
% Sarah scan were recorded in a different direction than in Montreal, so
% this code just swap the direction of the images for the processing.
% ======================================================================= %
if strcmp(bone_main_axis,'z')
elseif strcmp(bone_main_axis,'y')
    %% This is for having the bone main axis as the z axis
    % The canadian scans, and hence order was x, y, z, hence "1 2 3"
    % x, y, z
    % 1, 2, 3
    Image = int16(permute(Image,[2 3 1]));
elseif strcmp(bone_main_axis,'x')
    %% This is for having the other direction
    % z, x, y
    % 3, 1, 2
    Image = int16(permute(Image,[3 1 2]));
else
    error('The axis was not read properly');
end


end