function centroid_invasion_struct = spheroid_centroid_invasion(...
                                        directory,...
                                        pixel_files,...
                                        centroid_cell_id,...
                                        centroid_cell_type,...
                                        centroid_time,...
                                        centroid_xcom,...
                                        centroid_ycom,...
                                        centroid_zcom,...
                                        micron_scale,...
                                        twelve_hour_intervals,...
                                        centroid_invasion_struct,...
                                        struct_counter,...
                                        xdim,...
                                        ydim,...
                                        zdim...
                                        ) 


        spheroid_initial_conditions=[directory pixel_files(1).name];
        fid=fopen(spheroid_initial_conditions,'r');
        pixel_ic_data=textscan(fid,'%f %f %f %f %f %f','Delimiter','MultipleDelimsAsOne','HeaderLines',1);
        fclose(fid);
        spheroid_ic_x=pixel_ic_data{4};
        spheroid_ic_y=pixel_ic_data{5};
        spheroid_ic_z=pixel_ic_data{6};
        binary_spheroid_ic=logical(zeros(ydim,xdim,zdim));
        for I=1:length(spheroid_ic_x)
            binary_spheroid_ic(spheroid_ic_y(I)+1,spheroid_ic_x(I)+1,spheroid_ic_z(I)+1)=1;
        end


        
        
        
        

        time_index = ismember(centroid_time,twelve_hour_intervals);
        centroid_time=centroid_time(time_index);
        cell_id=centroid_cell_id(time_index);
        cell_type=centroid_cell_type(time_index);
        xcom=centroid_xcom(time_index);
        ycom=centroid_ycom(time_index);
        zcom=centroid_zcom(time_index);
        
        cell_type_index=find(cell_type==1);
        centroid_time_scc=centroid_time(cell_type_index);
        cell_id_scc=cell_id(cell_type_index);
        xcom_scc=xcom(cell_type_index);
        ycom_scc=ycom(cell_type_index);
        zcom_scc=zcom(cell_type_index);
        unique_scc=unique(cell_id_scc);
        
        scc_centroid_map=nan(max(unique_scc),length(twelve_hour_intervals));
        distance_transform_spheroid_ic=bwdist(binary_spheroid_ic);
        scc_centroid_location_index = sub2ind([ydim,xdim,zdim],round(ycom_scc(:))+1,round(xcom_scc(:))+1,round(zcom_scc(:))+1);
        distances = micron_scale * distance_transform_spheroid_ic(scc_centroid_location_index);

        for t=1:length(twelve_hour_intervals)
            t_index=find(centroid_time_scc==twelve_hour_intervals(t));
            input_map(:,:) = scc_centroid_map(:,t);
            input_map(cell_id_scc(t_index),1)=distances(t_index);
            scc_centroid_map(:,t) = input_map;
        end
        centroid_invasion_struct(struct_counter).invasion=scc_centroid_map;
end