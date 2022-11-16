function no_fractured_objects = fractured_object_counter(directory_output,no_combinations,no_repeats,quantification_day)

for combination_index=1:no_combinations
    mat_file=[directory_output 'neighbours_fractures_run' sprintf('%05lu',combination_index)];
    load(mat_file);
    for repeat_index=1:no_repeats
        removenan=isnan(fractured_object_struct(repeat_index,quantification_day).no_sccs);
        if(min(removenan)==0)
            no_invading_objects(repeat_index)=length(fractured_object_struct(repeat_index,quantification_day).no_sccs(~removenan));
        else
            no_invading_objects(repeat_index)=0;
        end;
        no_fractured_objects(combination_index,repeat_index)=no_invading_objects(repeat_index);
       
    end
        
 
end


   