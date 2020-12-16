
% for angle_diff = 1 : 180
%     
%     min_scl = 0;   % the min degree search range
%     gma = 0.015; % Motion explotion Average parameter
% 
%     score = (1*(((1-gma).^(angle_diff-1) + min_scl)/(1 + min_scl)));
%     
%     plot(angle_diff,score,'r.'); hold on;
% 
% end


function score = view_score(angle_diff)

    angle_diff = rad2deg(angle_diff);

    min_scl = 0;   % the min degree search range
    gma = 0.025; % Motion explotion Average parameter
    
    score = (1*(((1-gma).^(angle_diff-1) + min_scl)/(1 + min_scl)));
   
end