function cell_count = cell_counter(cell_count,centroid_time,centroid_cell_type,struct_counter,combination_index) 
index_caf=find(centroid_cell_type==2);
centroid_time(index_caf)=[];
       
unique_time=unique(centroid_time);
for t=1:length(unique_time)
    time_index=find(centroid_time==unique_time(t));
    cell_count(t,struct_counter,combination_index)=length(time_index);
end


        