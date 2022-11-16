function cell_growth_rates = cell_growth_rate(directory_output,no_combinations,no_repeats,quantification_day)

mat_file=fullfile(directory_output, 'cell_count.mat');
load(mat_file);
end_time = (quantification_day*2*60*24/40)+1;
t=0:end_time-1;
t=t*40/2;%Time in mins
cell_count_input=cell_count(1:end_time,:,:);

for combination_index=1:no_combinations
    for repeat_index=1:no_repeats
        single_cell_count(:)=cell_count_input(:,repeat_index,combination_index);
        p = polyfit(t,log(single_cell_count),1);
        cell_growth_rates(combination_index,repeat_index)=p(1);
    end
end
end



