function centroid_invasion_struct = organotypic_centroid_invasion(...
                                        centroid_cell_id,...
                                        centroid_cell_type,...
                                        centroid_time,...
                                        centroid_zcom,...
                                        invasion_boundary,...
                                        micron_scale,...
                                        twelve_hour_intervals,...
                                        centroid_invasion_struct,...
                                        struct_counter...
                                        ) 

if(isempty(centroid_zcom) == 0)
    centroid_zcom = - centroid_zcom + invasion_boundary; %Sets negative z to be no invasion
    scc_index=find(centroid_cell_type==1);
    time_scc=centroid_time(scc_index);
    cell_id_scc=centroid_cell_id(scc_index);
    zcom_scc=centroid_zcom(scc_index);
    unique_scc=unique(cell_id_scc);
    unique_scc_times=unique(time_scc);
    scc_centroid_map=nan(length(unique_scc),length(unique_scc_times));
    for t=1:length(unique_scc_times)
        time_index=find(time_scc==unique_scc_times(t));
        scc_centroid_map(1:length(time_index),t)=zcom_scc(time_index);
    end
    scc_centroid_map=micron_scale*scc_centroid_map;
    centroid_invasion_struct(struct_counter).invasion=scc_centroid_map(:,twelve_hour_intervals);
end
end