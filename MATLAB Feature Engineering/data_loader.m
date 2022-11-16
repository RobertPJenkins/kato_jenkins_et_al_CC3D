function [...
          directory,...
          centroid_time,...
          centroid_cell_id,...
          centroid_cell_type,...
          centroid_xcom,...
          centroid_ycom,...
          centroid_zcom,...
          neighbour_MC_step,...
          neighbour_cell_id,...
          neighbour_cell_type,...
          neighbour_cell_scc_neighbour,...
          neighbour_cell_caf_neighbour,...
          pixel_files...
          ]...
          = data_loader(directory_input,folders,correct_index,repeat_index)        

reordered_folder_index=correct_index(repeat_index);
folder_directory=folders(reordered_folder_index).name;
directory=[directory_input folder_directory filesep];

centroid_file_name=[directory '*CentroidData.txt'];
centroid_files=dir(centroid_file_name);
centroid_file=[directory centroid_files(1).name];
fid = fopen(centroid_file,'r');
centroid_data=textscan(fid,'%f %f %f %f %f %f','Delimiter','MultipleDelimsAsOne','HeaderLines',1);
fclose(fid);
centroid_time=centroid_data{1};
centroid_cell_id=centroid_data{2};
centroid_cell_type=centroid_data{3};
centroid_xcom=centroid_data{4};
centroid_ycom=centroid_data{5};
centroid_zcom=centroid_data{6};

neighbour_file_name=[directory '*NoNeighbours*'];
neighbour_files=dir(neighbour_file_name);
neighbour_file=[directory neighbour_files(1).name];
fid = fopen(neighbour_file,'r');
neighbour_data=textscan(fid,'%f %f %f %f %f','Delimiter','MultipleDelimsAsOne','HeaderLines',1);
fclose(fid);
neighbour_MC_step=neighbour_data{1};
neighbour_cell_id=neighbour_data{2};
neighbour_cell_type=neighbour_data{3};
neighbour_cell_scc_neighbour=neighbour_data{4};
neighbour_cell_caf_neighbour=neighbour_data{5};

pixel_file_names=[directory '*PixelData_MCS*'];
pixel_files=dir(pixel_file_names);

end