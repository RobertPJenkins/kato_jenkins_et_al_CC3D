function [...
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
                                            )

            pixel_file_t=[directory pixel_files(t+1).name];
            fid=fopen(pixel_file_t,'r');
            pixel_data_t=textscan(fid,'%f %f %f %f %f %f','Delimiter','MultipleDelimsAsOne','HeaderLines',1);
            fclose(fid);
            cell_id=pixel_data_t{2};
            cell_type=pixel_data_t{3};
            x=pixel_data_t{4};
            y=pixel_data_t{5};
            z=pixel_data_t{6};
            
            main_tissue_t(:,:,:) = main_tissue_all_t(:,:,:,t+1);
            
            MC_neighbour_t_index=find(neighbour_MC_step==time_vector(t+1));
            neighbour_cell_id_t=neighbour_cell_id(MC_neighbour_t_index);
            neighbour_cell_type_t=neighbour_cell_type(MC_neighbour_t_index);
            neighbour_cell_scc_neighbour_t=neighbour_cell_scc_neighbour(MC_neighbour_t_index);
            neighbour_cell_caf_neighbour_t=neighbour_cell_caf_neighbour(MC_neighbour_t_index);
            
            cells_from_pixels=zeros(ydim,xdim,zdim);
            cell_types_from_pixels=zeros(ydim,xdim,zdim);
            for cell_index=1:length(cell_id)
                cells_from_pixels(y(cell_index)+1,x(cell_index)+1,z(cell_index)+1)=cell_id(cell_index);
                cell_types_from_pixels(y(cell_index)+1,x(cell_index)+1,z(cell_index)+1)=cell_type(cell_index);
            end;
            main_tissue_cell_ids=cells_from_pixels.*double(main_tissue_t);
            main_tissue_cell_ids=main_tissue_cell_ids.*(~binary_main_tissue_ic);
            unique_main_tissue_cell_ids=unique(main_tissue_cell_ids(:));
            unique_main_tissue_cell_ids=unique_main_tissue_cell_ids(unique_main_tissue_cell_ids>0);
            
            index=find(centroid_time==time_vector(t+1));
            centroid_cell_id_t=centroid_cell_id(index);
            centroid_xcom_t=centroid_xcom(index);
            centroid_ycom_t=centroid_ycom(index);
            centroid_zcom_t=centroid_zcom(index);
            spheroid_centroid_cell_id_t=uint16(zeros(ydim,xdim,zdim));
            for cell_index=1:length(centroid_xcom_t)
                spheroid_centroid_cell_id_t(round(centroid_ycom_t(cell_index))+1,round(centroid_xcom_t(cell_index))+1,round(centroid_zcom_t(cell_index))+1)=centroid_cell_id_t(cell_index);
            end
            spheroid_centroid_cell_id_t=spheroid_centroid_cell_id_t.*uint16(~binary_main_tissue_ic);
            unique_centroid_cell_id_t=unique(spheroid_centroid_cell_id_t(:));
            unique_centroid_cell_id_t(unique_centroid_cell_id_t==0)=[];
            invasion_distance_t=bwdist(binary_main_tissue_ic);
            
            distance_index=find(spheroid_centroid_cell_id_t);
            spheroid_centroid_cell_id_t=spheroid_centroid_cell_id_t(distance_index);
            invasion_distance_t=invasion_distance_t(distance_index);
            
            unique_invasive_cell_ids_t=intersect(unique_centroid_cell_id_t,unique_main_tissue_cell_ids);





            fractured_cell_types_from_pixels=cell_types_from_pixels.*double(~main_tissue_t);
            fractured_cells_from_pixels=cells_from_pixels.*double(~main_tissue_t);

            fractured_objects=logical(fractured_cells_from_pixels);
            fractured_objects = bwareaopen(fractured_objects,150,18); %Fractured objects less than 150 voxels in size are then removed.
            fractured_object_labels=bwlabeln(fractured_objects,18);
            unique_fractured_objects=unique(fractured_object_labels);
            unique_fractured_objects=unique_fractured_objects(unique_fractured_objects>0);

            for object_index=1:length(unique_fractured_objects)
                label_index=find(fractured_object_labels==unique_fractured_objects(object_index));
                label_types=fractured_cell_types_from_pixels(label_index);
                unique_label_cell_types=unique(label_types);
                unique_label_cell_types=unique_label_cell_types(unique_label_cell_types>0);
                if(length(unique_label_cell_types)==1)
                   if( unique_label_cell_types==2)
                       fractured_object_labels(label_index)=0;
                   end
                end
            end;
            unique_fractured_objects=unique(fractured_object_labels);
            unique_fractured_objects=unique_fractured_objects(unique_fractured_objects>0);
end