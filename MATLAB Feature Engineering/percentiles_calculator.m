function [invasion_score_percentiles,maximum_invasion_percentiles,tapering_percentiles,mean_scc_neighbours_percentiles,fractured_object_percentiles,growth_rate_percentiles] = percentiles_calculator(invasion_score,maximum_invasion,tapering,mean_scc_neighbours,no_fractured_objects,cell_growth_rates,no_repeats,no_combinations)

for combination_index=1:no_combinations
    combination_invasion_score=zeros(no_repeats,1);
    combination_invasion_score(:,:)=invasion_score(combination_index,:);
    combination_maximum_invasion=zeros(no_repeats,1);
    combination_maximum_invasion(:,:)=maximum_invasion(combination_index,:);

    combination_tapering=zeros(no_repeats,1);
    combination_tapering(:,:)=tapering(combination_index,:);
    combination_mean_scc_neighbours=zeros(no_repeats,1);
    combination_mean_scc_neighbours(:,:)=mean_scc_neighbours(combination_index,:);
    combination_fractures=zeros(no_repeats,1);
    combination_fractures(:,:)=no_fractured_objects(combination_index,:);

    combination_count=zeros(no_repeats,1);
    combination_count(:,:)=cell_growth_rates(combination_index,:);

    for(i=1:19)
        invasion_score_percentiles(combination_index,i)=prctile(combination_invasion_score,5*i);
        maximum_invasion_percentiles(combination_index,i)=prctile(combination_maximum_invasion,5*i);
        tapering_percentiles(combination_index,i)=prctile(combination_tapering,5*i);
        mean_scc_neighbours_percentiles(combination_index,i)=prctile(combination_mean_scc_neighbours,5*i);
        fractured_object_percentiles(combination_index,i)=prctile(combination_fractures,5*i);
        growth_rate_percentiles(combination_index,i)=prctile(combination_count,5*i);
    end
end
end



