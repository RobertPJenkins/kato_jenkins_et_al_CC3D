function [invasion_stats_all_sccs,invasion_score,maximum_invasion] = invasion_stats_compiler(directory_output,no_combinations,no_repeats,quantification_day)

time_point=2*quantification_day;
for combination_index=1:no_combinations
    mat_file=[directory_output 'invasion_centroids_run' sprintf('%05lu',combination_index)];
    load(mat_file);
    for repeat_index=1:no_repeats
        single_scc_invasion=centroid_invasion_struct(repeat_index).invasion;
        [~,T]=size(single_scc_invasion);
        for t=1:T%Time of event
            single_scc_invasion_t = single_scc_invasion(:,t);
            single_scc_invasion_t=single_scc_invasion_t(~isnan(single_scc_invasion_t));
            single_scc_invasion_t=single_scc_invasion_t(single_scc_invasion_t>0);
            %Number of invading objects.
            invasion_stats_single_scc(1,t,repeat_index)=length(single_scc_invasion_t);
            if(length(single_scc_invasion_t)>0)
                invasion_stats_single_scc(2,t,repeat_index)=mean(single_scc_invasion_t(:));
                invasion_stats_single_scc(3,t,repeat_index)=max(single_scc_invasion_t(:));
            else %All set to zero if no invasion.
                invasion_stats_single_scc(2,t,repeat_index)=0;
                invasion_stats_single_scc(3,t,repeat_index)=0;
            end
        end
    end
invasion_stats_all_sccs(:,:,:,combination_index)=invasion_stats_single_scc;
metric=zeros(no_repeats,1);
metric(:)=invasion_stats_all_sccs(1,time_point+1,:,combination_index).*invasion_stats_all_sccs(2,time_point+1,:,combination_index);
invasion_score(combination_index,:) = metric;
metric(:)=invasion_stats_all_sccs(3,time_point+1,:,combination_index);
maximum_invasion(combination_index,:) = metric;
  
end



end

