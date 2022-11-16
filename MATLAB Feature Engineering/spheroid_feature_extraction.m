function spheroid_feature_extraction(directory_input,directory_output,no_combinations,no_repeats,micron_scale,xdim,ydim,zdim,quantification_day)

mkdir(directory_output)
folder_names=[directory_input '*run*'];
folders=dir(folder_names);
isub = [folders(:).isdir];
index=find(isub==0);
folders(index)=[];
twelve_hour_intervals=0:1440:14400;
time_vector=[0,2880,5760,8640,11520,14400];

correct_index = folder_reordering(no_combinations,no_repeats,folders);
cell_count=[];
for combination_index=1:no_combinations

    struct_counter=0;
    centroid_invasion_struct(no_repeats)=struct('invasion', []);
    main_tissue_struct(no_repeats) = struct('tissue',[]);
    neighbours_struct(no_repeats,5)=struct('invading_cell_neighbours', []);
    fractured_object_struct(no_repeats,5)=struct('no_sccs',[],'no_cafs',[],'min_object_depth',[],'max_object_depth',[]);

    scc_struct(no_repeats)=struct();
    
    for repeat_index=no_repeats*(combination_index-1)+1:no_repeats*combination_index
        struct_counter=struct_counter+1;
        combination_index, struct_counter
        [...
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
          = data_loader(directory_input,folders,correct_index,repeat_index); 
        
        cell_count = cell_counter(cell_count,centroid_time,centroid_cell_type,struct_counter,combination_index);

        centroid_invasion_struct = spheroid_centroid_invasion(...
                                                        directory,...
                                                        pixel_files,...
                                                        centroid_cell_id,...
                                                        centroid_cell_type,...
                                                        centroid_time,...
                                                        centroid_xcom,...
                                                        centroid_ycom,...
                                                        centroid_zcom,...
                                                        micron_scale,...
                                                        twelve_hour_intervals,...
                                                        centroid_invasion_struct,...
                                                        struct_counter,...
                                                        xdim,...
                                                        ydim,...
                                                        zdim...
                                                        );


        [main_tissue_all_t, binary_main_tissue_ic] = ...
                                    spheroid_main_tissue_extractor(...
                                        directory,...
                                        pixel_files,...
                                        xdim,...
                                        ydim,...
                                        zdim,...
                                        time_vector,...
                                        neighbour_MC_step,...
                                        neighbour_cell_id,...
                                        neighbour_cell_scc_neighbour,...
                                        neighbour_cell_caf_neighbour,...
                                        main_tissue_struct,...
                                        struct_counter...
                                        );

         for t=quantification_day
    [...
        spheroid_centroid_cell_id_t,...
        unique_invasive_cell_ids_t,...
        neighbour_cell_id_t,...
        neighbour_cell_scc_neighbour_t,...
        neighbour_cell_caf_neighbour_t,...
        neighbour_cell_type_t,...
        invasion_distance_t,...
        unique_fractured_objects,...
        fractured_object_labels,...
        fractured_cells_from_pixels,...
        fractured_cell_types_from_pixels...
    ] = spheroid_neighbour_fracture_data_loader(...
                                        directory,...
                                        pixel_files,...
                                        main_tissue_all_t,...
                                        binary_main_tissue_ic,...
                                        struct_counter,...
                                        centroid_time,...
                                        centroid_cell_id,...
                                        centroid_xcom,...
                                        centroid_ycom,...
                                        centroid_zcom,...
                                        neighbour_MC_step,...
                                        neighbour_cell_id,...
                                        neighbour_cell_type,...
                                        neighbour_cell_scc_neighbour,...
                                        neighbour_cell_caf_neighbour,...
                                        t,...
                                        time_vector,...
                                        xdim,...
                                        ydim,...
                                        zdim...
                                            );


    neighbours_struct = spheroid_neighbour_structure_compiler(...
                                        spheroid_centroid_cell_id_t,...
                                        unique_invasive_cell_ids_t,...
                                        invasion_distance_t,...
                                        neighbour_cell_id_t,...
                                        neighbour_cell_type_t,...
                                        neighbour_cell_scc_neighbour_t,...
                                        neighbour_cell_caf_neighbour_t,...
                                        centroid_cell_id,...
                                        neighbours_struct,...
                                        struct_counter,...
                                        t...
                                        );   

fractured_object_struct = ...
    spheroid_fractured_object_compiler(...
                             unique_fractured_objects,...
                             fractured_object_labels,...
                             fractured_cells_from_pixels,...
                             fractured_cell_types_from_pixels,...
                             ydim,...
                             xdim,...
                             zdim,...
                             micron_scale,...
                             fractured_object_struct,...
                             struct_counter,...
                             t...
                             );
        end
        mat_file_save_name=fullfile(directory_output, ['main_tissue_run_' sprintf('%05lu',correct_index(repeat_index)) '.mat']);
        save(mat_file_save_name,'main_tissue_all_t','combination_index','repeat_index','correct_index');
    end

    mat_file_save_name=fullfile(directory_output, ['neighbours_fractures_run' sprintf('%05lu',combination_index) '.mat']);
    save(mat_file_save_name,'neighbours_struct','fractured_object_struct','combination_index');
    
    mat_file_save_name=fullfile(directory_output, ['invasion_centroids_run' sprintf('%05lu',combination_index) '.mat']);
    save(mat_file_save_name,'centroid_invasion_struct','combination_index');






end

mat_file_save_name=fullfile(directory_output,'cell_count.mat');
save(mat_file_save_name,'cell_count');


[invasion_stats_all_sccs,invasion_score,maximum_invasion] = invasion_stats_compiler(directory_output,no_combinations,no_repeats,quantification_day);
[tapering,mean_scc_neighbours] = spheroid_strand_statistics(directory_output,no_combinations,no_repeats,micron_scale,zdim,quantification_day);
no_fractured_objects = fractured_object_counter(directory_output,no_combinations,no_repeats,quantification_day);
% cell_growth_rates = cell_growth_rate(directory_output,no_combinations,no_repeats,quantification_day);
cell_growth_rates = cell_growth_rate(directory_output,no_combinations,no_repeats,quantification_day);
[invasion_score_percentiles,maximum_invasion_percentiles,tapering_percentiles,mean_scc_neighbours_percentiles,fractured_object_percentiles,growth_rate_percentiles] = percentiles_calculator(invasion_score,maximum_invasion,tapering,mean_scc_neighbours,no_fractured_objects,cell_growth_rates,no_repeats,no_combinations);
mat_file_save_name=fullfile(directory_output, 'compiled_data.mat');
save(mat_file_save_name,'invasion_stats_all_sccs','invasion_score','maximum_invasion','tapering','mean_scc_neighbours','no_fractured_objects','cell_growth_rates','invasion_score_percentiles','maximum_invasion_percentiles','tapering_percentiles','mean_scc_neighbours_percentiles','fractured_object_percentiles','growth_rate_percentiles');

end