function fractured_object_metrics = spheroid_fractured_object_quantifier(unique_fractured_objects,fractured_object_labels,fractured_cells_from_pixels,fractured_cell_types_from_pixels,ydim,xdim,zdim)

for fractured_object_index=1:length(unique_fractured_objects)
   index=find(fractured_object_labels==unique_fractured_objects(fractured_object_index));
   [y_coords,x_coords,z_coords]=ind2sub([ydim,xdim,zdim],index);
   
   y_coords=abs(y_coords- round(ydim/2));
   x_coords=abs(x_coords-round(xdim/2));
   z_coords=abs(z_coords-round(zdim/2));
   
   single_object_cells=fractured_cells_from_pixels(index);
   unique_single_object_cells=unique(single_object_cells);
   for cell_index=1:length(unique_single_object_cells)
       index=find(fractured_cells_from_pixels==unique_single_object_cells(cell_index));
       cell_type_counter(cell_index)=unique(fractured_cell_types_from_pixels(index));
   end
   fractured_object_metrics(fractured_object_index,1)=length(find(cell_type_counter==1));
   fractured_object_metrics(fractured_object_index,2)=length(find(cell_type_counter==2));
   fractured_object_metrics(fractured_object_index,3)=min([min(x_coords),min(y_coords),min(z_coords)]);
   fractured_object_metrics(fractured_object_index,4)=max([max(x_coords),max(y_coords),max(z_coords)]);
end

end


