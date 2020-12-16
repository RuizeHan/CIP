function [match_scores,data] = VecMatching(data)
    
    vec_top = data.vec_top;  % normalizaing£¿
    vec_ego = data.vec_ego;
    [img_h,img_w,img_z] = size(data.img_ego);
   
    factor = 5;

    search_num = length(vec_top);
    match_scores = zeros(1,search_num);
    match_idx = cell(1,search_num);
    match_way = zeros(1,search_num);
    
    for i =  1 : search_num
        
        vec1 = vec_top{i};  % normalizaing£¿
        vec2 = vec_ego;
        len2 = length(vec2(1,:));
        
        
        if size(vec1,2) < 2
%             match_score(i) = Inf;
%             continue;
              match_way(i) = -1;
              match_scores(i) = factor^len2;
        else
            
        len1 = length(vec1(1,:));
        len2 = length(vec2(1,:));
        %vec1(1,:) = normalizing(vec1(1,:),0,1);
        vec1(1,:) = vec1(1,:);
        vec1(2,:) = normalizing(vec1(2,:),0,1);
        vec2(1,:) = vec2(1,:)/(img_w/2);
        vec2(2,:) = normalizing(-vec2(2,:),0,1);
        
%         filepath = 'D:\Research\multi_camera\vis_vec_top\';
%         filename = [filepath num2str(i),'.jpg'];
%         h = figure(6);
%         plot(vec1(1,:),vec1(2,:),'ro','Linewidth',1.5);
%         title(['The search angle is : ',num2str(2*i),'¡ã']);
%         axis([-1 1 0 max(vec1(2,:))])
%         saveas(h,filename);


        if  len1 == len2       
            match_scores(i) = sum((vec1(:) - vec2(:)).^2)/len1; % x_1==0  
            match_way(i) = 0;
        elseif len1 < len2
            [match_scores(i),match_idx{i}] = get_match_score(vec1,vec2,factor);    
            match_way(i) = -1;
        else
            [match_scores(i),match_idx{i}] = get_match_score(vec2,vec1,factor); 
            match_way(i) =1;
        end
        
        end
    end
    
    data.match_idx = match_idx;
    data.match_way = match_way;
    
end
