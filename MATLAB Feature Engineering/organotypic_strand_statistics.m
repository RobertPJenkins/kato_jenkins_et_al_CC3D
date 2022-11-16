function [tapering,mean_scc_neighbours] = organotypic_strand_statistics(directory_output,no_combinations,no_repeats,micron_scale,zdim,tissue_boundary,quantification_day)

for combination_index=1:no_combinations

    mat_file=[directory_output 'neighbours_fractures_run' sprintf('%05lu',combination_index)];
    load(mat_file);

    for repeat_index=1:no_repeats
        cell_type=neighbours_struct(repeat_index,quantification_day).invading_cell_neighbours(:,2);
        scc_neigh=neighbours_struct(repeat_index,quantification_day).invading_cell_neighbours(:,3);
        depth=neighbours_struct(repeat_index,quantification_day).invading_cell_neighbours(:,5);
        
        %Consider only depths greater than the tissue boundary
        index=find(depth>zdim-tissue_boundary);
        scc_neigh = scc_neigh(index);
        cell_type=cell_type(index);
        depth=depth(index) - (zdim-tissue_boundary);
        depth = micron_scale * depth;

        index=intersect(find(cell_type==1),find(scc_neigh>0));
        scc_neigh=scc_neigh(index);
        depth=depth(index);

        p=polyfit(depth,scc_neigh,1);
        tapering(combination_index,repeat_index)=p(1);
        mean_scc_neighbours(combination_index,repeat_index)=mean(scc_neigh);
    end
end
