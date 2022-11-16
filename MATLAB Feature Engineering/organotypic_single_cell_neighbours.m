function invading_cell_neighbours = organotypic_single_cell_neighbours(...
    unique_possible_neighbours,...
    neighbour_cell_id_t,...
    neighbour_cell_type_t,...
    neighbour_cell_scc_neighbour_t,...
    neighbour_cell_caf_neighbour_t,...
    centroid_cell_id,...
    centroid_zcom,...
    centroid_time,...
    input_time)


counter=0;
for neighbour_index=1:length(unique_possible_neighbours)
    unique_neighbour_index=find(neighbour_cell_id_t==unique_possible_neighbours(neighbour_index));
    
    if(length(unique_neighbour_index)>0)
        no_scc_neighbours=neighbour_cell_scc_neighbour_t(unique_neighbour_index);
        no_caf_neighbours=neighbour_cell_caf_neighbour_t(unique_neighbour_index);
        
        if(no_scc_neighbours>0)||(no_caf_neighbours>0)
            counter=counter+1;
            invading_cell_neighbours(counter,1)=neighbour_cell_id_t(unique_neighbour_index);
            invading_cell_neighbours(counter,2)=neighbour_cell_type_t(unique_neighbour_index);
            invading_cell_neighbours(counter,3)=no_scc_neighbours;
            invading_cell_neighbours(counter,4)=no_caf_neighbours;
        
            time_index=find(centroid_time==input_time);
            neighbour_centroids=centroid_cell_id(time_index);
            neighbour_centroid_zcom=centroid_zcom(time_index);
        
            neighbour_index=find(neighbour_centroids==unique_possible_neighbours(neighbour_index));
            depth=neighbour_centroid_zcom(neighbour_index);
            invading_cell_neighbours(counter,5)=depth;
        end
    else
        invading_cell_neighbours(counter,1)=NaN;
        invading_cell_neighbours(counter,2)=NaN;
        invading_cell_neighbours(counter,3)=NaN;
        invading_cell_neighbours(counter,4)=NaN;
        invading_cell_neighbours(counter,5)=NaN;
    end

    
end
if(length(unique_possible_neighbours)==0)
    invading_cell_neighbours(1,1)=NaN;
    invading_cell_neighbours(1,2)=NaN;
    invading_cell_neighbours(1,3)=NaN;
    invading_cell_neighbours(1,4)=NaN;
    invading_cell_neighbours(1,5)=NaN;
end
end
