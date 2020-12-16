function focus_level = single_focus(params,candidate,subjects,can_i)
    
    focus_level = cell(1,3);
    angle_diff = zeros(length(subjects),1);
    view_type = zeros(length(subjects),1);
    count_focus = 0;
    focus_id = zeros(1,length(subjects)-1);
%     figure(8) 
%     imshow(data.img_top)    
%     hold on;
    for view_i = setdiff(1:length(subjects),can_i)   %  All the other persons observe the 'candidate'
        
        connect_diff_x = candidate.x - subjects{view_i}.x;
        connect_diff_y = candidate.y - subjects{view_i}.y;
        con = [connect_diff_x,connect_diff_y];
        
        a = con;
        view_start = [subjects{view_i}.x,subjects{view_i}.y]; 
        view_angle = subjects{view_i}.view_angle;
        view_confi = subjects{view_i}.view_confi;
        
        if view_confi == 0
           view_type(view_i) = 0;
        elseif view_confi >=1
           view_type(view_i) = 0.5;
        else
           view_type(view_i) = 1;
        end

        L = 50; 
        % view_end = [view_start(1) + L * cos(view_angle) ,  view_start(2) + L * sin(view_angle)];
        view_end = [view_start(1) + L * view_angle(1) ,  view_start(2) - L * view_angle(2)];
        b = view_end-view_start;
        
%         figure(8)
%         plot(candidate.x,candidate.y,'rh');
%         plot(view_start(1),view_start(2),'ro');
%         text(view_start(1)+10,view_start(2)+10,num2str(view_i),'Color','yellow','FontSize',15)
%         plot([view_start(1),view_end(1)],[view_start(2),view_end(2)],'g');
            
        
        angle_differ = real(acos(dot(a,b)/(norm(a)*norm(b))));%
        angle_differ(isnan(angle_differ)==1) = 0;
        angle_diff(view_i) =  angle_differ;
        % angle_diff(view_i) = subspace(a',b');  % The view_angle difference between the other viewers and one candidate
        angle_diff_deg = rad2deg(angle_diff);
       
        if angle_diff(view_i) < params.focus_threshold && any(view_angle)
            count_focus = count_focus + 1;
            focus_id(count_focus) = view_i;
        end
        
    end
    
    view_scores = view_type .* view_score(angle_diff); % view_score
    
    focus_level{1} = sum(view_scores)/sum(view_type~=0);  % The average focus_level
    
    focus_level{2} = count_focus;     % The totle focusing persons
    focus_level{3} = focus_id(1:count_focus);

end