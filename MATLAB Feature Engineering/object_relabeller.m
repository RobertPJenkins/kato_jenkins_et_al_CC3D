function [unique_object_labels,object_label] = object_relabeller(object_label)
unique_object_labels=unique(object_label);
index=find(unique_object_labels==0);
unique_object_labels(index)=[];

counter=0;
for object_index=1:length(unique_object_labels)
    counter=counter+1;
    index=find(object_label==unique_object_labels(object_index));
    object_label(index)=counter;
end
unique_object_labels=unique(object_label);
index=find(unique_object_labels==0);
unique_object_labels(index)=[];