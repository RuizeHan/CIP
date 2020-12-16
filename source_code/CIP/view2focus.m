clear
clf
clc

X = [10,10,50,50,65,80,75];
Y= [10,50,50,10,65,55,40];

View_angle_1 = [45,-0,-135,135,-135,-135,135]; %  Right : 0  Left : 180
View_angle_2 = [45,-0,-135,135,-135,-135,135]; %  Right : 0  Left : 180
View_angle_3 = [45,-0,-135,135,-135,-135,135]; %  Right : 0  Left : 180


View_angle_fusion = [45,-0,-135,135,-135,-135,135]; %  Right : 0  Left : 180

View_angle = View_angle_fusion * 2 * pi/360;

subjects = cell(length(X),1);

params.focus_threshold = pi/6;

figure(1)
hold on;
axis([0 100 0 100]);

for can_i = 1 : length(X)
    
    candi.x = X(can_i);
    candi.y = Y(can_i);
    candi.view_angle = View_angle(can_i);
    subjects{can_i} = candi;
    
    view_start = [subjects{can_i}.x,subjects{can_i}.y]; 
    view_angle = subjects{can_i}.view_angle;

    L = 10;  % The visual length of the view_angle

    view_end = [view_start(1) + L * cos(view_angle) ,  view_start(2) + L * sin(view_angle)];
    plot(view_start(1),view_start(2),'rs');
    text(view_start(1)+1,view_start(2)+1,num2str(can_i))
    plot([view_start(1),view_end(1)],[view_start(2),view_end(2)]);
    
end

focus_level = cell(length(subjects),3);

for can_i = 1:length(subjects)
   
    candidate =  subjects{can_i};
    
    focus_level_i = single_focus(params,candidate,subjects,can_i);
    
    focus_level(can_i,:) = focus_level_i;
    
end


overlap_t_r = zeros(length(subjects_t),length(subjects_t));
focus_smooth_t_r = zeros(length(subjects_t),length(subjects_t));

for can_t_i = 1:length(subjects_t)  %  All candidates from frame_t
    
    candidate_t_i =  subjects_t{can_t_i};
    focus_level_t_i = single_focus(params,candidate_t_i,subjects_t,can_t_i); 
        
    for can_r_i = 1:length(subjects_r)   %  All candidates from frame_r
 
        %focus_level_t(can_t_i,:) = focus_level_t_i;

        candidate_r_i =  subjects{can_r_i};
        focus_level_r_i = single_focus(params,candidate_r_i,subjects_r,can_r_i);  
        %focus_level_r(can_i,:) = focus_level_r_i;
        
        focus_smooth_t_r(can_t_i,can_r_i) = focus_level_r_i - focus_level_t_i;

        overlap_t_r(can_t_i,can_r_i) = get_overlap(candidate_t_i,candidate_r_i);
 
    end
    
end







