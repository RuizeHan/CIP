function res_gaze = generate_gaze(res,data,params)
    
    vis = params.vis_gaze;
    img_top = data{1}.img_top;
    [h,w,~] = size(img_top);
    I = img_top;
    if vis == 1
       figure(7)
       imshow(I);
       set (gcf,'Position',[300,100,600,400])
       hold on;
       text(100,100,num2str(params.framenum),'color','y','fontsize',25);
    end
    gopro_num = size(res,1);
    color_space = ['c','m','y'];
    
    gaze_all = cell(gopro_num,size(data{1}.objs_top,1));
    
    objs_top = data{1}.objs_top;
    obj_x = objs_top(:,1)+1/2*objs_top(:,3);
    obj_y = objs_top(:,2)+1/2*objs_top(:,4);
    
    if vis == 1
        for sub_i = 1 : length(obj_x)
            plot(obj_x(sub_i),obj_y(sub_i),'.g','Linewidth',5);
            text(obj_x(sub_i)-30,obj_y(sub_i)-30,num2str(sub_i),'color','g','fontsize',15);
        end
    end
    
    for gopro_i = 1:gopro_num
        res_gopro_i = res{gopro_i};
        data_gopro_i = data{gopro_i};
        cam_top = data_gopro_i.cam_top;
        cam = [cam_top(1)+1/2*cam_top(3),cam_top(2)+1/2*cam_top(4)];
            
        if vis == 1
            plot(cam(1),cam(2),'pr','Linewidth',5);
  
            vertical_ang = res_gopro_i.search_angle;
            right_ang = vertical_ang + pi/4;
            left_ang = vertical_ang + pi/4*3;
            V = [cos(vertical_ang),-sin(vertical_ang)]; 
            LV  = get_line_equation(cam(1),cam(2),V(1)+cam(1),V(2)+cam(2));
            if LV(2) == 0
                plot([-LV(3)/LV(1),-LV(3)/LV(1)],[0,h],color_space(gopro_i),'LineWidth',1.5);               %%%%%%plot angle
            else
                plot([0,w],[-LV(3)/LV(2),(-LV(3)-LV(1)*w)/LV(2)],color_space(gopro_i),'LineWidth',1.5);      %%%%%%plot angle
            end
            
            corner_x = [w, 0, 0, w];
            corner_y = [0, 0, h, h];
            cor_dx = corner_x - cam(1);
            cor_dy = corner_y - cam(2);

            % get the angle of the line contected the camera with each corner 
            num_cor =  length(corner_x);
            angle_cor = zeros(1,num_cor);
            for i = 1 : num_cor
                vec_x = [1,0];
                vec_cor = [cor_dx(i),cor_dy(i)];
                angle_cor(i) = acos(dot(vec_x,vec_cor)/(norm(vec_x)*norm(vec_cor)));
                if cor_dy(i) > 0
                    angle_cor(i) = 2 * pi - angle_cor(i);
                end
            end
            for judge_angle_ori  = [left_ang,right_ang]
               judge_angle = mod(judge_angle_ori,2*pi);
               if judge_angle >= angle_cor(1) && judge_angle < angle_cor(2)
                    L1 = [0,1,cam(2)];
                elseif judge_angle >= angle_cor(2) && judge_angle < angle_cor(3)
                    L1 = [1,0,cam(1)];
                elseif judge_angle >= angle_cor(3) && judge_angle < angle_cor(4)
                    L1 = [0,1,cam(2)-h];
                elseif judge_angle <= angle_cor(1) || judge_angle > angle_cor(4)
                    L1 = [1,0,cam(1)-w];
               end
                L2  = get_line_equation(0,0,cos(judge_angle),-sin(judge_angle));
                [X,Y] = get_cross_point(L1,L2);
                if vis ==1
                    plot([ceil(X+cam(1)),cam(1)],[ceil(Y+cam(2)),cam(2)],color_space(gopro_i),'LineWidth',1.5);
                end
            end
        end

        idx_ego_ord = res_gopro_i.ego.ord;
        idx_ego_eff = res_gopro_i.ego.eff;

        idx_top_eff = res_gopro_i.top.eff;
        idx_top_ord = res_gopro_i.top.ord;

        objs_ego = data_gopro_i.objs_ego;
        objs_top = data_gopro_i.objs_top;

        search_angle = res_gopro_i.search_angle;
        search_angle = search_angle * 180/pi; 

        gaze_line_h = 50;
        
        for i = 1:length(idx_top_eff)
            
        idx_top_obj = idx_top_eff(idx_top_ord(i));

        %--------------------------------------------------------------
        
        center_x = objs_top(idx_top_eff(idx_top_ord(i)), 1) + 0.5 * objs_top(idx_top_eff(idx_top_ord(i)), 3);
        center_y = objs_top(idx_top_eff(idx_top_ord(i)), 2) + 0.5 * objs_top(idx_top_eff(idx_top_ord(i)), 4);
        
        connect_diff_x = center_x - cam(1);
        connect_diff_y = center_y - cam(2);
        con = [connect_diff_x,connect_diff_y];
        
        view_angle = search_angle;
        view_start = [cam(1),cam(2)]; 
        view_end = [view_start(1) + cosd(view_angle), view_start(2) - sind(view_angle)];

        a = con;
        b = view_end-view_start;
        
        angle_diff = real(acos(dot(a,b)/(norm(a)*norm(b))));%
        
        hor_view_dir = objs_ego(idx_ego_ord(idx_ego_eff(i)), 5);
        
        %hor_view_dir_trans = mod( - hor_view_dir - 90, 360);  % From Down : 0 to Right : 0 
       
        gaze_ang = - hor_view_dir + search_angle + rad2deg(angle_diff) + 180;
           
        %----------------------------------------------------------------
         
        gaze_angle = mod(deg2rad(gaze_ang), 2*pi);   % From Down : 0 to Right : 0 
     
        angle_confidence = objs_ego(idx_ego_ord(idx_ego_eff(i)), 6);

           objs_top(idx_top_eff(idx_top_ord(i)), 5) = gaze_angle;

           center_x = objs_top(idx_top_eff(idx_top_ord(i)), 1) + 0.5 * objs_top(idx_top_eff(idx_top_ord(i)), 3);
           center_y = objs_top(idx_top_eff(idx_top_ord(i)), 2) + 0.5 * objs_top(idx_top_eff(idx_top_ord(i)), 4);
           
%          right 0  anticlockwise
%          if gaze_angle > -90 && gaze_angle < 90
           gaze_x = center_x + gaze_line_h * cos(gaze_angle) * angle_confidence;
           gaze_y = center_y - gaze_line_h * sin(gaze_angle) * angle_confidence;
           
           if vis == 1
                plot([center_x, gaze_x], [center_y, gaze_y],'Linewidth',1.5,'color',color_space(gopro_i));
           end
           
           gaze_all{gopro_i,idx_top_obj} = [gaze_angle,angle_confidence];
           
        end
        
        data_gopro_i.objs_top = objs_top;
        
    end
    
    %  compute the view_angle of the camera wearer
    for gopro_i = 1:gopro_num  
        
        data_gopro_i = data{gopro_i};
        cam_idx = data_gopro_i.cam_id;
        
        res_gopro_i = res{gopro_i};
        search_angle = res_gopro_i.search_angle;
        view_angle = search_angle + pi/2;
        angle_confidence = 1;
        
        for i = 1:gopro_num  
            gaze_all{i,cam_idx} = [view_angle,angle_confidence];
        end
        
    end

    res_gaze = gaze_all;
    
    for fus_type = params.fus_type
        gaze_fusion = angle_fusion(gaze_all,fus_type);
    end
    
    subjects = cell(size(objs_top,1),1);
    
    for sub_i = 1 : size(objs_top,1)
   
       gaze_fuse  = gaze_fusion{sub_i};
       center_x = objs_top(sub_i, 1) + 0.5 * objs_top(sub_i, 3);
       center_y = objs_top(sub_i, 2) + 0.5 * objs_top(sub_i, 4);
       gaze_x = center_x + gaze_line_h * gaze_fuse(1);
       gaze_y = center_y - gaze_line_h * gaze_fuse(2);
       
       if vis == 1
           figure(7)
           plot([center_x, gaze_x], [center_y, gaze_y],'Linewidth',1.5,'color','r'); 
       end
       
       candi.x = center_x;
       candi.y = center_y;
       candi.objs_top = objs_top(sub_i,1:4);
       
       candi.view_angle = gaze_fuse;
       
       candi.view_confi = sqrt(gaze_fuse(1)^2+gaze_fuse(2)^2); 
       
       subjects{sub_i} = candi;

    end

    res_gaze = subjects;
end