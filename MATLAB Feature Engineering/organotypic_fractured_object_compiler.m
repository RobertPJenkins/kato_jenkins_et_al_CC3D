function fractured_object_struct = ...
    organotypic_fractured_object_compiler(...
                             fractured_object_struct,...
                             binary_all_cells_from_pixels,...
                             cells_from_pixels,...
                             cell_types_from_pixels,...
                             index_all_tissue,...
                             xdim,...
                             ydim,...
                             zdim,...
                             micron_scale,...
                             struct_counter,...
                             t...
                             )

binary_all_cells_from_pixels(index_all_tissue)=0;
fractured_cell_label = bwlabeln(binary_all_cells_from_pixels,18);
[unique_fractured_cells,fractured_cell_label] = join_periodic_objects(fractured_cell_label,ydim,xdim);
%remove objects only with cafs
for object_index=1:length(unique_fractured_cells)
    index=find(fractured_cell_label==unique_fractured_cells(object_index));
    object_type=cell_types_from_pixels(index);
    unique_object_type=unique(object_type(:));
    if(length(unique_object_type)==1)
        if(unique_object_type==2)
            fractured_cell_label(index)=0;
            binary_all_cells_from_pixels(index)=0;
        end;
    end;
end;
unique_fractured_cells=unique(fractured_cell_label);
index_all_tissue=find(unique_fractured_cells==0);
unique_fractured_cells(index_all_tissue)=[];
%Fractured objects less than 150 voxels in size are then removed.
for object_index=1:length(unique_fractured_cells)
    index=find(fractured_cell_label==unique_fractured_cells(object_index));
    if(length(index)<150)
        fractured_cell_label(index)=0;
    end
end;
[unique_fractured_cells,fractured_cell_label] = object_relabeller(fractured_cell_label);

%Here we collect fractured object info
clear fractured_object_metrics
if(length(unique_fractured_cells)>0)
    fractured_object_metrics = organotypic_fractured_object_quantifier(fractured_cell_label,cells_from_pixels,cell_types_from_pixels,ydim,xdim,zdim);
    fractured_object_struct(struct_counter,t).no_sccs=fractured_object_metrics(:,1);
    fractured_object_struct(struct_counter,t).no_cafs=fractured_object_metrics(:,2);
    fractured_object_struct(struct_counter,t).min_object_depth=micron_scale*fractured_object_metrics(:,3);
    fractured_object_struct(struct_counter,t).max_object_depth=micron_scale*fractured_object_metrics(:,4);
else
    fractured_object_struct(struct_counter,t).no_sccs=NaN;
    fractured_object_struct(struct_counter,t).no_cafs=NaN;
    fractured_object_struct(struct_counter,t).min_object_depth=NaN;%In microns
    fractured_object_struct(struct_counter,t).max_object_depth=NaN;%In microns
end

end