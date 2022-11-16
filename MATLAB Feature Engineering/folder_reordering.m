function correct_index = folder_reordering(no_combinations,no_repeats,folders)

for i=1:no_combinations*no_repeats
     folder_names{i}=folders(i).name;
end;
current_index=1:no_combinations*no_repeats;
counter=0;
for i=1:no_combinations
    for j=1:no_repeats
        counter=counter+1;
        correct_index(counter)=current_index((j-1)*no_combinations+i);
    end
end