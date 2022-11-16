function fractured_object_metrics = organotypic_fractured_object_quantifier(fractured_cell_label,cells_from_pixels,cell_types_from_pixels,ydim,xdim,zdim)

unique_fractured_cells=unique(fractured_cell_label);
index=find(unique_fractured_cells==0);
unique_fractured_cells(index)=[];

    for object_index=1:length(unique_fractured_cells)
       index=find(fractured_cell_label==unique_fractured_cells(object_index));
       [~,~,z_label_pixels]=ind2sub([ydim,xdim,zdim],index);
       fractured_object=cells_from_pixels(index);
       unique_object_cells=unique(fractured_object);
       cell_type_counter = object_counter(cells_from_pixels,cell_types_from_pixels,unique_object_cells);
       fractured_object_metrics(object_index,1)=length(find(cell_type_counter==1));%No. SCCs
       fractured_object_metrics(object_index,2)=length(find(cell_type_counter==2));%No. CAFs
       fractured_object_metrics(object_index,3)=min(z_label_pixels);%Depth min
       fractured_object_metrics(object_index,4)=max(z_label_pixels);%Depth max
    end
end


function cell_type_counter = object_counter(cells_from_pixels,cell_types_from_pixels,unique_object_cells)
for J=1:length(unique_object_cells)
   index=find(cells_from_pixels==unique_object_cells(J));
   cell_type_counter(J)=unique(cell_types_from_pixels(index));
end
end
