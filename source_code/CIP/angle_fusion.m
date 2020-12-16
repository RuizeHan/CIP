function gaze_fusion = angle_fusion(gaze_all,fus_type)

gaze_fusion = cell(size(gaze_all,2),1);

    for sub_i = 1 : size(gaze_all,2)
        
        gaze_fus = [0 0];
    
        gaze = cell2mat(gaze_all(:,sub_i));
            
         if ~isempty(gaze)
                        
            if fus_type == 1
                
                max_con_gaze = gaze((gaze(:,2) == max(gaze(:,2))),:);

                if size(max_con_gaze,1) > 1
                    gaze_ang = max_con_gaze(1,1);
                    gaze_cof = max_con_gaze(1,2);
                else
                    gaze_ang = max_con_gaze(1);
                    gaze_cof = max_con_gaze(2);
                end

                gaze_vector = gaze_cof * [cos(gaze_ang),sin(gaze_ang)];

    %             figure(11)
    %             axis([-1 1 -1 1]);
    %             plot([0 gaze_vector(1)],[0 gaze_vector(2)]);
    %             hold on;

                gaze_fus = gaze_vector;
            
            elseif fus_type ==2   % vector add
                
                viewers = size(gaze_all,1);
                for view_i = 1 : viewers

                    gaze = gaze_all{view_i,sub_i};

                    if ~isempty(gaze)

                    gaze_ang = gaze(1);
                    gaze_cof = 1/viewers * gaze(2);

                    gaze_vector = gaze_cof * [cos(gaze_ang),sin(gaze_ang)];

                    gaze_fus = gaze_fus + gaze_vector;

                    end

                end
                   
            else   % single view gaze
                switch fus_type
                    case -1
                        if ~isempty(gaze(1,:))
                        max_con_gaze = gaze(1,:);
                        end
                        
                    case -2
                        if ~isempty(gaze(2,:))
                        max_con_gaze = gaze(2,:);
                        end
                        
                    case -3
                        if ~isempty(gaze(3,:))
                        max_con_gaze = gaze(3,:);
                        end
                    
                end
                gaze_ang = max_con_gaze(1);
                gaze_cof = max_con_gaze(2);
                gaze_vector = gaze_cof * [cos(gaze_ang),sin(gaze_ang)];
                gaze_fus = gaze_vector; 
            end
            
        end
       
%         plot([0 gaze_fus(1)],[0 gaze_fus(2)]);
%         hold off;
        
    gaze_fusion{sub_i} = gaze_fus;

    end


end