function save_visualization_full(which_eval, day, str_path, Cortical_both, Image_t0, Image_tX, ...
    Result_volumes_inclsurf_labels, Result_volumes_inclsurf, constant_bone, formed_bone, resorbed_bone, ...
    Endosteal_both, Periosteal_both, dist_formed, loc_max_formed, dist_resorbed, loc_max_resorbed)
% ======================================================================= %
%% 2022-03-10
% save_visualization_full(which_eval, day, str_path, Cortical_both, Image_t0, Image_tX, ...
%     Result_volumes_inclsurf_labels, Result_volumes_inclsurf, constant_bone, formed_bone, resorbed_bone, ...
%     Endosteal_both, Periosteal_both, dist_formed, loc_max_formed, dist_resorbed, loc_max_resorbed)
% This will create all of the outputs and save them as 3D tiff stacks.
% ======================================================================= %
Day_cell = num2str(day);
% saving the visualization one by one
boneboth = strcat(which_eval,'Both', Day_cell);
imwrite3D_stacked(Cortical_both, str_path, boneboth);


it0 = strcat('Image_t0_', Day_cell);
imwrite3D_stacked(Image_t0, str_path, it0);


itX = strcat('Image_tX_', Day_cell);
imwrite3D_stacked(Image_tX, str_path, itX);


RVISL = strcat('Result_volumes_inclsurf_labels_', Day_cell);
imwrite3D_stacked(Result_volumes_inclsurf_labels, str_path, RVISL);


RVIS = strcat('Result_volumes_inclsurf_', Day_cell);
imwrite3D_stacked(Result_volumes_inclsurf, str_path, RVIS);

constantbone = strcat('constant_bone_', Day_cell);
imwrite3D_stacked(constant_bone, str_path, constantbone);


formedbone = strcat('formed_bone_', Day_cell);
imwrite3D_stacked(formed_bone, str_path, formedbone);


resorbedbone = strcat('resorbed_bone_', Day_cell);
imwrite3D_stacked(resorbed_bone, str_path, resorbedbone);


switch which_eval
    case 'cort'
        Endostealboth = strcat('Endosteal_both_', Day_cell);
        imwrite3D_stacked(Endosteal_both, str_path, Endostealboth);
        
        Periostealboth = strcat('Periosteal_both_', Day_cell);
        imwrite3D_stacked(Periosteal_both, str_path, Periostealboth);
    case 'trab'
end


distformed = strcat('dist_formed_', Day_cell);
imwrite3D_stacked(dist_formed, str_path, distformed);


locmaxformed = strcat('loc_max_formed_', Day_cell);
imwrite3D_stacked(loc_max_formed, str_path, locmaxformed);


distres = strcat('dist_resorbed_', Day_cell);
imwrite3D_stacked(dist_resorbed, str_path, distres);


locmaxres = strcat('loc_max_resorbed_', Day_cell);
imwrite3D_stacked(loc_max_resorbed, str_path, locmaxres);

end