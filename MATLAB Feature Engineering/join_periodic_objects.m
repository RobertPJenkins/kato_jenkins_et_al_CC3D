function [unique_object_labels,output_label] = join_periodic_objects(input_label,ydim,xdim)  

periodic_objects = padarray(input_label,[ydim,xdim,0],'circular');

periodic_object_labels=bwlabeln(logical(periodic_objects),18);
unique_periodic_objects=unique(periodic_object_labels(:));
unique_periodic_objects(unique_periodic_objects==0)=[];
output_label = input_label;
for periodic_object_index=1:length(unique_periodic_objects);
    index=find(periodic_object_labels==unique_periodic_objects(periodic_object_index));
    overlap=unique(periodic_objects(index));
    if(length(overlap)>1)
        for J=2:length(overlap)
            index=find(output_label==overlap(J));
            output_label(index)=overlap(1);
        end;
    end;
end
unique_object_labels=unique(output_label);
index=find(unique_object_labels==0);
unique_object_labels(index)=[];