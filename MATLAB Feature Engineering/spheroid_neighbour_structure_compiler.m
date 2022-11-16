function neighbours_struct = spheroid_neighbour_structure_compiler(...
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
    )

counter=0;
for cell_index=1:length(unique_invasive_cell_ids_t)
    cell_id_index=find(neighbour_cell_id_t==unique_invasive_cell_ids_t(cell_index));
    
    if(length(cell_id_index)>0)
        no_scc_neighbours=neighbour_cell_scc_neighbour_t(cell_id_index);
        no_caf_neighbours=neighbour_cell_caf_neighbour_t(cell_id_index);
        
        if(no_scc_neighbours>0)||(no_caf_neighbours>0)
            counter=counter+1;
            invading_cell_neighbours(counter,1)=neighbour_cell_id_t(cell_id_index);
            invading_cell_neighbours(counter,2)=neighbour_cell_type_t(cell_id_index);
            invading_cell_neighbours(counter,3)=no_scc_neighbours;
            invading_cell_neighbours(counter,4)=no_caf_neighbours;
            
            depth_index=find(spheroid_centroid_cell_id_t==unique_invasive_cell_ids_t(cell_index));
            invading_cell_neighbours(counter,5)=invasion_distance_t(depth_index);
        else
            invading_cell_neighbours(counter,1)=NaN;
            invading_cell_neighbours(counter,2)=NaN;
            invading_cell_neighbours(counter,3)=NaN;
            invading_cell_neighbours(counter,4)=NaN;
            invading_cell_neighbours(counter,5)=NaN;
        end
    end
end
if(length(unique_invasive_cell_ids_t)==0)
    invading_cell_neighbours(1,1)=NaN;
    invading_cell_neighbours(1,2)=NaN;
    invading_cell_neighbours(1,3)=NaN;
    invading_cell_neighbours(1,4)=NaN;
    invading_cell_neighbours(1,5)=NaN;
end
neighbours_struct(struct_counter,t).invading_cell_neighbours=invading_cell_neighbours;
end





