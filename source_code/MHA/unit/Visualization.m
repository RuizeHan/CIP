function res = Visualization(objs_top,cam_top,objs_ego,data,params,flag)
    
    vis = params.vis;

    img_top = data.img_top;
    img_ego = data.img_ego;
    idx_eff_top = data.index_top{data.max_idx};
    idx_eff_ego  = data.index_ego;
        
    if flag ==1
        if data.match_way(data.max_idx) == 1  % length of vec_top < length of vec_ego
            idx_eff_top  = idx_eff_top(data.match_idx{data.max_idx}); 
        elseif data.match_way(data.max_idx) == -1 % length of vec_top > length of vec_ego
            idx_eff_ego  = data.match_idx{data.max_idx};     
        end
    else
        if isempty(idx_eff_top)
            idx_eff_top = 1;
        else
            idx_eff_top  = idx_eff_top(data.match_idx_top{data.max_idx}); 
            idx_eff_ego  = idx_eff_ego(data.match_idx_ego{data.max_idx}); 
        end
    end
    
    
    obj_x = objs_ego(:,1);
    obj_y = objs_ego(:,2);
    obj_w = objs_ego(:,3);
    obj_h = objs_ego(:,4);

    I = img_ego;
    obj_num = length(obj_x);
    
    cet_x = obj_x + 0.5 * obj_w;
    cet_y = obj_y + 0.5 * obj_h;
         
    [cet_x_ord,idx_ego_ord] = sort(cet_x);
    cet_y_ord = cet_y(idx_ego_ord);
    
    cet_x_eff = cet_x_ord(idx_eff_ego);
    cet_y_eff = cet_y_ord(idx_eff_ego);
    
    res.ego.ord = idx_ego_ord;
    res.ego.eff = idx_eff_ego;
    
    if vis == 1
        figure(3)
        imshow(I);
        set (gcf,'Position',[300,100,600,400])

        hold on;
        plot(cet_x,cet_y,'ro','Linewidth',1.5);


        for i = 1:obj_num   
            rectangle('Position',[obj_x(i),obj_y(i),obj_w(i),obj_h(i)],'LineWidth',2,'EdgeColor','b');        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         gaze = textread('F:\MHA\gaze_dataset\gaze_direction.txt');
%         for i = 1:size(gaze, 2)/5
%             plot([gaze(1 + 4*(i - 1)), gaze(3 + 4*(i - 1))], [gaze(2 + 4*(i - 1)), gaze(4 + 4*(i - 1))],'Linewidth',1.5,'color','r');
%         end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
% % % % %  hor - gt
%         for i = 1:size(data.objs_ego_gt,1)
%             rectangle('Position',[data.objs_ego_gt(i,1),data.objs_ego_gt(i,2),data.objs_ego_gt(i,3),data.objs_ego_gt(i,4)],'LineWidth',2,'EdgeColor','g');        
%         end
% 
%         for i = 1 : size(data.lab.ego,1) 
%             text(data.objs_ego_gt(i,1),data.objs_ego_gt(i,2),num2str(data.lab.ego(i)),'color','g','fontsize',15);
%         end 

        for i = 1 : length(cet_x_eff) 
            text(cet_x_eff(i),cet_y_eff(i),num2str(i),'color','r','fontsize',15);
        end
    
    end
    %---------------------------------UAV-----------------------------------

    [h,w,z] = size(img_top);
    cam = [cam_top(1)+1/2*cam_top(3),cam_top(2)+1/2*cam_top(4)];
    obj_x = objs_top(:,1)+1/2*objs_top(:,3);
    obj_y = objs_top(:,2)+1/2*objs_top(:,4);
    I = img_top;

    obj_num = length(obj_x);
    
    if vis == 1
        figure(1)
        imshow(I);
        set (gcf,'Position',[300,600,600,400])  
        hold on;
        plot(cam(1),cam(2),'pr','Linewidth',5);
    %    plot(data.objs_top_gt(:,1)+1/2*data.objs_top_gt(:,3),data.objs_top_gt(:,2)+1/2*data.objs_top_gt(:,4),'og','Linewidth',1.5)
    %     for i = 1:size(data.objs_top_gt,1)
    %         rectangle('Position',[data.objs_top_gt(i,1),data.objs_top_gt(i,2),data.objs_top_gt(i,3),data.objs_top_gt(i,4)],'LineWidth',2,'EdgeColor','r');
    %     end

        for i = 1:obj_num   
            rectangle('Position',[objs_top(i,1),objs_top(i,2),objs_top(i,3),objs_top(i,4)],'LineWidth',2,'EdgeColor','b');        
        end
        text(100,100,num2str(params.framenum),'color','y','fontsize',25);
    
    end
    
    angle = data.angle_objs{data.max_idx};
    search_angle = data.current_view_angle;
    vertical_ang = search_angle;
    left_ang = vertical_ang + pi/4*3;
    
    angle_eff = angle(idx_eff_top);
    if left_ang < 2 * pi
        angle_eff(angle_eff > 2*pi) = angle_eff(angle_eff > 2*pi) - 2*pi;
    end

    [~,idx_ord] = sort(angle_eff,'descend');
    obj_x_eff = obj_x(idx_eff_top);
    obj_y_eff = obj_y(idx_eff_top);
    obj_x_ord = obj_x_eff(idx_ord);
    obj_y_ord = obj_y_eff(idx_ord);
    
    res.top.eff = idx_eff_top;
    res.top.ord = idx_ord;
    res.search_angle = search_angle;
    
    if vis == 1
        for i = 1:length(obj_x_eff)
            text(obj_x_ord(i),obj_y_ord(i),num2str(i),'color','r','fontsize',15);
        end 
        
%%%%%%%%%%% top - gt
%         for i = 1 : size(data.lab.top,1)
%             text(data.objs_top_gt(i,1)+50,data.objs_top_gt(i,2)+50,num2str(data.lab.top(i)),'color','g','fontsize',15);
%         end 

        obj_dx = obj_x - cam(1);
        obj_dy = obj_y - cam(2);
        angle = zeros(1,obj_num);

        % get the angle of each object
        for i = 1 : obj_num  
            vec_x = [1,0];
            vec_obj = [obj_dx(i),obj_dy(i)];
            angle(i) = acos(dot(vec_x,vec_obj)/(norm(vec_x)*norm(vec_obj)));
            if obj_dy(i) > 0
                angle(i) = 2 * pi - angle(i);
            end
        end
        angle_deg = angle * 180/pi;

        search_num =params.search_num;
        %   vec_top = cell(search_num,1);
        for k = data.max_idx
%           search_angle = pi/(search_num/2) * (k);            
            search_angle = data.current_view_angle;
            vertical_ang = search_angle;
            right_ang = vertical_ang + pi/4;
            left_ang = vertical_ang + pi/4*3;

            % get the distance of each object in specfic view 
            flag = zeros(1,obj_num);
            dis = zeros(1,obj_num);
            foot = zeros(2,obj_num);

            for i = 1 : obj_num  

                if left_ang > 2 * pi  % !todo : when the leftangle greater than 2*pi      
                    ang_gre = find(right_ang < 1/2*pi & right_ang > 0);
                    angle(ang_gre) = angle(ang_gre) + 2 * pi;        
                end

               P = [obj_dx(i),obj_dy(i)];
               O = [0,0];    
               V = [cos(vertical_ang),-sin(vertical_ang)];  
               OV = [0,0,cos(vertical_ang),-sin(vertical_ang)];  

                if angle(i) > right_ang && angle(i) < left_ang     % !todo : when the leftangle greater than 2*pi
                   flag(i) = 1;   
                   proj_point = get_foot_point(P,OV);
                   foot(1,i) = proj_point(1);
                   foot(2,i) = proj_point(2);
                   if vis ==1
                        %plot(proj_point(1)+cam(1),proj_point(2)+cam(2),'rs');
                   end

                   %  Option 1 : the distance of the objects P to the vertical line LV
                   % dis(i) = abs(det([V-O;P-O]))/norm(V-O);
                   %  Option 2 : the distance of the objects P to the camera point O
                   dis(i) = norm(P-O);

                end    
            end

            % get the equation of vertical line in the original axis
            LV  = get_line_equation(cam(1),cam(2),V(1)+cam(1),V(2)+cam(2));

            if vis ==1
                if LV(2) == 0
                    plot([-LV(3)/LV(1),-LV(3)/LV(1)],[0,h],'c','LineWidth',1.5);               %%%%%%plot angle
                else
                    plot([0,w],[-LV(3)/LV(2),(-LV(3)-LV(1)*w)/LV(2)],'c','LineWidth',1.5);      %%%%%%plot angle
                end
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

          % get the cross point of the side
          for judge_angle_ori  = [left_ang,right_ang]
               judge_angle = mod(judge_angle_ori,2*pi);
                % L = [cos(left_ang),-sin(left_ang)];
                % R = [cos(right_ang),-sin(right_ang)]; 
                % judge which edge to cross  
                edge_index = 0;
                if judge_angle >= angle_cor(1) && judge_angle < angle_cor(2)
                    edge_index = 1;
                    L1 = [0,1,cam(2)];
                elseif judge_angle >= angle_cor(2) && judge_angle < angle_cor(3)
                    edge_index = 2;
                    L1 = [1,0,cam(1)];
                elseif judge_angle >= angle_cor(3) && judge_angle < angle_cor(4)
                    edge_index = 3; 
                    L1 = [0,1,cam(2)-h];
                elseif judge_angle <= angle_cor(1) || judge_angle > angle_cor(4)
                    edge_index = 4;
                    L1 = [1,0,cam(1)-w];
                end

                L2  = get_line_equation(0,0,cos(judge_angle),-sin(judge_angle));
                [X,Y] = get_cross_point(L1,L2);
                if vis ==1
                    plot([ceil(X+cam(1)),cam(1)],[ceil(Y+cam(2)),cam(2)],'b','LineWidth',1.5);
                end
          end
        end   
    end

end
