function [data] = UAV2vec_fast(img,obj,cam,params,data)

    [h,w,z] = size(img);
    cam = [cam(1)+1/2*cam(3),cam(2)+1/2*cam(4)];
    obj_x = obj(:,1)+1/2*obj(:,3);
    obj_y = obj(:,2)+1/2*obj(:,4);
    I =img;

    vis = 0;
    obj_num = length(obj_x);
    
    if vis == 1
        figure(1)
        set (gcf,'Position',[300,600,600,400])
        imshow(I);
        hold on;
        plot(cam(1),cam(2),'pg','Linewidth',5);
        plot(obj_x(:),obj_y(:),'or','Linewidth',1.5);
    end

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
 
    
    angle_deg = angle * 180/pi; %弧度制，转角度制乘180/pi

%     search_num = params.search_num;
    for k =  1 : 2*obj_num
        
        temp = mod(k,2);
        temp(temp==0)=3;
        
        search_angle = angle(ceil(k/2)) - temp*pi/4;
        
        if search_angle<=0
            search_angle = search_angle+2*pi;
        end

        vertical_ang = search_angle;
        right_ang = vertical_ang + pi/4;
        left_ang = vertical_ang + pi/4*3;

        % get the distance of each object in specfic view 
        flag = zeros(1,obj_num);
        dis = zeros(1,obj_num);
        foot = zeros(2,obj_num);
        
        

        for i = 1 : obj_num  

            if left_ang > 2 * pi  % !todo : when the leftangle greater than 2*pi      
                ang_gre = find(angle < 3*pi/4 & angle > 0);
                angle(ang_gre) = angle(ang_gre) + 2 * pi;        
            end

           P = [obj_dx(i),obj_dy(i)];
           O = [0,0];    
           V = [cos(vertical_ang),-sin(vertical_ang)];  
           OV = [0,0,cos(vertical_ang),-sin(vertical_ang)];  

            if angle(i) >= right_ang && angle(i) <= left_ang     % !todo : when the leftangle greater than 2*pi
               flag(i) = 1;
               proj_point = get_foot_point(P,OV);
               foot(1,i) = proj_point(1);
               foot(2,i) = proj_point(2);
               %   if vis ==1
               %        %plot(proj_point(1)+cam(1),proj_point(2)+cam(2),'rs');
               %   nd
               
               %  Option 1 : the distance of the objects P to the vertical line LV
               dis(i) = abs(det([V-O;P-O]))/norm(V-O);
               %  Option 2 : the distance of the objects P to the camera point O
%              dis(i) = norm(P-O);               
            end    
        end
        
        % Occlusion
        if params.occ == 1
        for m = 1 : obj_num - 1
            for n = m + 1 : obj_num
                if  abs(angle_deg(m) - angle_deg(n)) < 2 % obj is occluded when the angle is approciate with another     
                    if dis(m) >= dis(n)
                        flag(m) = 0;
                    else
                        flag(n) = 0; 
                    end
                end
            end
        end
        end
        
        angle_diff = angle - vertical_ang;
        
        theta_ind = vertical_ang * 180/pi;
        theta_ind = round(theta_ind/(params.interval));
        if theta_ind==0
            theta_ind=round(360/params.interval);
        end
        
        [data.vec_top{theta_ind},data.index_top{theta_ind}] = get_distribution_vector(flag,angle_diff,foot,dis);
        data.angle = angle;

    end
    

end
