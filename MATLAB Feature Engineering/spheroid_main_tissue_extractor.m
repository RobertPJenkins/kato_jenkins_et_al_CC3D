function [main_tissue_all_t, binary_main_tissue_ic] = spheroid_main_tissue_extractor(...
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
                                                )
spheroid_ic_file=[directory pixel_files(1).name];
fid=fopen(spheroid_ic_file,'r');
pixel_ic_data=textscan(fid,'%f %f %f %f %f %f','Delimiter','MultipleDelimsAsOne','HeaderLines',1);
fclose(fid);
spheroid_ic_cell_id=pixel_ic_data{2};
spheroid_ic_x=pixel_ic_data{4};
spheroid_ic_y=pixel_ic_data{5};
spheroid_ic_z=pixel_ic_data{6};
cell_id_pixels_ic=zeros(ydim,xdim,zdim);
binary_spheroid_ic=logical(zeros(ydim,xdim,zdim));
for pixel_index=1:length(spheroid_ic_x)
    cell_id_pixels_ic(spheroid_ic_y(pixel_index)+1,spheroid_ic_x(pixel_index)+1,spheroid_ic_z(pixel_index)+1)=spheroid_ic_cell_id(pixel_index);
    binary_spheroid_ic(spheroid_ic_y(pixel_index)+1,spheroid_ic_x(pixel_index)+1,spheroid_ic_z(pixel_index)+1)=1;
end

binary_spheroid_ic=bwareaopen(binary_spheroid_ic,50,26);
object_labels_ic = bwlabeln(binary_spheroid_ic,18);
unique_object_labels_ic=unique(object_labels_ic);
index=find(unique_object_labels_ic==0);
unique_object_labels_ic(index)=[];
for object_index=1:length(unique_object_labels_ic)
    object_size(object_index)=length(find(object_labels_ic==unique_object_labels_ic(object_index)));
end

main_tissue_cell_id=(zeros(ydim,xdim,zdim));
index1=find(object_size==max(object_size));
index2=find(object_labels_ic==unique_object_labels_ic(index1));
main_tissue_cell_id(index2)=cell_id_pixels_ic(index2);
unique_main_tissue_cell_ids=unique(main_tissue_cell_id(:));
unique_main_tissue_cell_ids=unique_main_tissue_cell_ids(unique_main_tissue_cell_ids>0);
index=find(neighbour_MC_step==0);
neighbours_ic_cell_id=neighbour_cell_id(index);
neighbours_ic_scc_neighbour=neighbour_cell_scc_neighbour(index);
neighbours_ic_caf_neighbour=neighbour_cell_caf_neighbour(index);
neighbours_ic_total_neighbour=neighbours_ic_caf_neighbour+neighbours_ic_scc_neighbour;

for I=1:length(unique_main_tissue_cell_ids)
    main_tissue_neighbour_index=find(neighbours_ic_cell_id==unique_main_tissue_cell_ids(I));
    if(neighbours_ic_total_neighbour(main_tissue_neighbour_index)==0)
        cell_remove_index=find(main_tissue_cell_id==unique_main_tissue_cell_ids(I));
        main_tissue_cell_id(cell_remove_index)=0;
    end;
end

binary_main_tissue_ic=logical(main_tissue_cell_id);
main_tissue_all_t(:,:,:,1)=binary_main_tissue_ic;



for t=1:length(time_vector)-1
    spheroid_t_file=[directory pixel_files(t+1).name];
    fid=fopen(spheroid_t_file,'r');
    pixel_t_data=textscan(fid,'%f %f %f %f %f %f','Delimiter','MultipleDelimsAsOne','HeaderLines',1);
    fclose(fid);
    cell_id=pixel_t_data{2};
    spheroid_t_x=pixel_t_data{4};
    spheroid_t_y=pixel_t_data{5};
    spheroid_t_z=pixel_t_data{6};
    cell_id_pixels_t=zeros(ydim,xdim,zdim);       
    binary_spheroid_t=logical(zeros(ydim,xdim,zdim));
    for pixel_index=1:length(cell_id)
        cell_id_pixels_t(spheroid_t_y(pixel_index)+1,spheroid_t_x(pixel_index)+1,spheroid_t_z(pixel_index)+1)=cell_id(pixel_index);
        binary_spheroid_t(spheroid_t_y(pixel_index)+1,spheroid_t_x(pixel_index)+1,spheroid_t_z(pixel_index)+1)=1;
    end
    
    object_labels_t = bwlabeln(binary_spheroid_t,18);
    main_tissue_object_labels_t=double(binary_main_tissue_ic).*object_labels_t;
    unique_main_tissue_object_labels_t=unique(main_tissue_object_labels_t(:));
    unique_main_tissue_object_labels_t=unique_main_tissue_object_labels_t(unique_main_tissue_object_labels_t>0);
    for object_index=1:length(unique_main_tissue_object_labels_t)
        index=find(object_labels_t==unique_main_tissue_object_labels_t(object_index));
        main_tissue_object_labels_t(index)=unique_main_tissue_object_labels_t(object_index);
    end

    cell_id_pixels_t=double(logical(main_tissue_object_labels_t)).*cell_id_pixels_t;
    unique_cell_id_pixels_t=unique(cell_id_pixels_t);
    unique_cell_id_pixels_t=unique_cell_id_pixels_t(unique_cell_id_pixels_t>0);
    
    index=find(neighbour_MC_step==time_vector(t+1));
    neighbours_t_cell_id=neighbour_cell_id(index);
    neighbours_t_scc_neighbour=neighbour_cell_scc_neighbour(index);
    neighbours_t_caf_neighbour=neighbour_cell_caf_neighbour(index);
    neighbours_t_total_neighbour=neighbours_t_caf_neighbour+neighbours_t_scc_neighbour;
    
    for cell_index=1:length(unique_cell_id_pixels_t)
        main_tissue_neighbour_index=find(neighbours_t_cell_id==unique_cell_id_pixels_t(cell_index));
        if(neighbours_t_total_neighbour(main_tissue_neighbour_index)==0)
            cell_remove_index=find(main_tissue_cell_id==unique_cell_id_pixels_t(cell_index));
            cell_id_pixels_t(cell_remove_index)=0;
        end;
    end

    binary_spheroid_t=logical(cell_id_pixels_t); 
    binary_spheroid_t=bwareaopen(binary_spheroid_t,50,26);
    main_tissue_all_t(:,:,:,t+1)=binary_spheroid_t;

end
end














