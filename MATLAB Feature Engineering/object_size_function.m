function object_sizes = object_size_function(object_label)
unique_object_labels=unique(object_label);
index=find(unique_object_labels==0);
unique_object_labels(index)=[];
object_sizes=[];
for object_index=1:length(unique_object_labels)
    object_sizes(object_index)=length(find(object_label==unique_object_labels(object_index)));
end