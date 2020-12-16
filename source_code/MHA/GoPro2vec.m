function [vec_ego,index_ego] = GoPro2vec(img,obj)
    
    vis = 0;
    [img_h,img_w,img_z] = size(img);
    obj_x = obj(:,1);
    obj_y = obj(:,2);
    obj_w = obj(:,3);
    obj_h = obj(:,4);

    % w = 300;
    % h = 200;
    % obj_x = [20,60,150,200,250];
    % obj_y = [80,50,10,80,50];
    % obj_w = [30,25,20,40,30];
    % obj_h = [100,75,45,120,100];

    I = img;
    obj_num = length(obj_x);
    index_ego = 1:obj_num;
    if vis == 1
        figure(3)
        set (gcf,'Position',[300,100,600,400])
        imshow(I);
        hold on;

        for i = 1:obj_num
            rectangle('Position',[obj_x(i),obj_y(i),obj_w(i),obj_h(i)],'LineWidth',2,'EdgeColor','r');
        end
    end
    cet_x = obj_x + 0.5 * obj_w;
    cet_y = obj_y + 0.5 * obj_h;
    
    if vis == 1
        plot(cet_x,cet_y,'ro','Linewidth',1.5);
    end

    vec_ego = get_distribution_vector_land(cet_x,obj_h,img_w);

%   figure(4)
%   plot(vec_ego(1,:),vec_ego(2,:),'ro');
%   axis([-inf inf 0 max(vec_ego(2,:))])

end
