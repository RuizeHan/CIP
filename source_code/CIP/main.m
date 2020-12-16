% Code for Complementary-View Co-Interest Person Detection

% By Ruize Han; Jiewen Zhao; Yiyang Gan
% End of the year 2019   
% Published in ACM Multimedia (ACM MM) 2020

close all
clear
clc
warning off 
addpath('../CIP');
addpath('../MHA/');
addpath('../MHA/unit');

params.new_gaze_res = 1;
params.new_CAP_res = 1;
params.vis_CAP_res = 1;
params.vis = 0;
params.vis_gaze = 0;
params.focus_threshold = deg2rad(45);
params.motion_dir = 0;   % Use the motion direction for look-at direction
params.fus_type = 1;

% % % % -------------------------------parameter settings------------------
params = params_setting(params);
scenenum = [1,2,3];

all_seq = configSeqs;
l = 1;
for i = 1 : length(scenenum)
    s_num = 3*(scenenum(i)-1)+1;
    all_scene(l:l+2) = all_seq(s_num:s_num+2);
    l = l + 3;
end

count=0;
times = 0;

for scenenum_i = 1 : 3 : length(all_scene)
 
    filePath = 'G:\Raiser\Research\CIP-MM2020\CIP_Dataset\';
    params.data_directory = 'G:\Raiser\Research\CIP-MM2020\CIP_public\CIP\';
    
    params.scene_i = all_scene{scenenum_i};
    img_path = [filePath,params.scene_i.name];
    Files = dir(fullfile(img_path,'/hor/*.jpg'));

    LengthFiles = length(Files);

    data.last_angles = zeros(params.num_gopro,1);
    all_sub_gaze = cell(LengthFiles,1);
    all_res_gaze = cell(LengthFiles,1);
    clip_length = 50;
    clip_num = ceil(LengthFiles/clip_length);
    video_res = ones(LengthFiles,1)*-1;
    bestC = zeros(1,clip_num);
    candidateOfClip = zeros(2,clip_num);

if params.new_gaze_res == 1    

% err_frame = zeros(LengthFiles,1);

for clip_i = 1 : clip_num

    frame_range = (clip_i-1) * clip_length + 1 : min((clip_i-1) * clip_length + clip_length,LengthFiles); 

    for frame_clip_i = 1 : length(frame_range)

        frame = frame_range(frame_clip_i);
        params.frame_clip_i = frame_clip_i;
        params.frame = frame;

        if params.temp_angle == 1
            min_scl = params.min_scl_degree;
            N_deg = floor(360*(((1-params.gma).^(frame-1) + min_scl)/(1 + min_scl)));
            params.search_num = ceil(N_deg/params.interval);
        else
            params.search_num = ceil(360/params.interval);
        end

        res_gaze = cell(params.num_gopro,1);
        data_gaze = cell(params.num_gopro,1);

        for gopro_i = 1 : params.num_gopro
            
            data.gopro_id = gopro_i;
            scenenum_gopro_i = scenenum_i + gopro_i - 1;
            if params.search_gopro == 2
                gopro_trk = load([filePath,'Tra_Gopro/',all_scene{scenenum_gopro_i}.name,'.mat']);
                gopro_trk = gopro_trk.results.res;
            end
            
            img_path = [filePath,all_scene{scenenum_gopro_i}.name];

            params.framenum = num2str(frame,'%04d');
            params.scenenum = all_scene{scenenum_gopro_i}.name;

            disp(['scene:',params.scenenum,' clip:',num2str(clip_i),' frame:',num2str(frame_clip_i),' allframe:',params.framenum]);

            img_ego_name = [params.framenum,'.jpg'];
            img_top_name = [params.framenum,'.jpg'];

            data.img_ego = imread([img_path,'/hor/',img_ego_name]);
            data.img_top = imread([img_path(1:end-2),'T','/top/',img_top_name]);

            [data.objs_top,data.objs_ego] = frame2obj(filePath,params);
 
            data.vec_top = cell(params.search_num,1);
            data.angle_objs = cell(params.search_num,1);
            data.index_top = cell(params.search_num,1);
            data.cam_top = gopro_trk(frame,:);
            
            cam_diff = (data.objs_top(:,1:2) - data.cam_top(1:2));
                        
            [min_dis,data.cam_id] = min((sum(abs(cam_diff),2)));
            data.objs_top_ngopro = data.objs_top(setdiff(1:size(data.objs_top,1),data.cam_id),:); % remove the hor camera in top view

            [data.vec_ego,data.index_ego] = GoPro2vec(data.img_ego,data.objs_ego);    

            data = UAV2vec(data.img_top,data.objs_top_ngopro,data.cam_top,params,data);     

            [match_scores,data] = VecMatching_Ransec_DTW(data,params);

            [max_sco,data.max_idx] = min(match_scores);

            % visualization for top_idx & ego_idx
            scale_y = data.scales_y(data.max_idx);

            % data.current_view_angle = mod(data.last_angles(gopro_i) + pi/(params.search_num/2) * (data.max_idx),2*pi);
            data.current_view_angle = mod(data.last_angles(gopro_i) - deg2rad(params.search_num/2) + deg2rad(data.max_idx),2*pi);
            
            res = Visualization(data.objs_top_ngopro,data.cam_top,data.objs_ego,data,params,2);

            data.last_angles(gopro_i) = data.current_view_angle;
            
            % call back
            hor_ref_idx = 1:size(data.objs_ego,1);
            hor_idx_ord = hor_ref_idx(res.ego.ord);
            res.hor_idx_corr = hor_idx_ord(res.ego.eff);
            
            top_ref_idx = setdiff(1:size(data.objs_top,1),data.cam_id);
            top_idx_eff = top_ref_idx(res.top.eff);
            res.top_idx_corr = top_idx_eff(res.top.ord);
            
%-----------------show image---------------------------
%             figure(5)
%             imshow(data.img_ego);
%             for i = 1:size(data.objs_ego,1)
%                 hold on;
%                 obj = data.objs_ego(i,:);
%                 cet_x = obj(1);
%                 cet_y = obj(2);
%                 plot(cet_x,cet_y,'ro','Linewidth',1.5);
%                 text(cet_x,cet_y,num2str(i),'color','r','fontsize',15);
%             end
%             
%             figure(6)
%             imshow(data.img_top);
%             for i = 1:size(data.objs_top,1)
%                 hold on;
%                 obj = data.objs_top(i,:);
%                 cet_x = obj(1);
%                 cet_y = obj(2);
%                 plot(cet_x,cet_y,'ro','Linewidth',1.5);
%                 text(cet_x,cet_y,num2str(i),'color','r','fontsize',15);
%             end
%-----------------show image---------------------------

            res_gaze{gopro_i} = res;
            data_gaze{gopro_i} = data;

        end
        
        all_res_gaze{frame} = res_gaze;  % association results of all subjects
        subjects_gaze = generate_gaze(res_gaze, data_gaze, params);  % gaze results of all subjects
        all_sub_gaze{frame} = subjects_gaze; % gaze results of all frames

    end   

end
    %% %% Save the subject and corresponding gaze results
    save(fullfile(params.data_directory,'Data',['res_' params.scene_i.scene '.mat']),'all_res_gaze','-v7.3');
    save(fullfile(params.data_directory,'Data',['gaze_' params.scene_i.scene '.mat']),'all_sub_gaze','-v7.3');
    % save(fullfile(params.data_directory,'Gaze',['err_' params.scene_i.scene '.mat']),'err_frame','-v7.3');

else
    load(fullfile(params.data_directory,'Data',['gaze_' params.scene_i.scene '.mat']));
    if params.motion_dir == 1
       all_sub_gaze = motion_direction(all_sub_gaze,params.scene_i);     % motion results of all subjects (ablation study)      
    end

end

if params.new_CAP_res == 1

    for clip_i = 1 : clip_num

    frame_range = (clip_i-1) * clip_length + 1 : min((clip_i-1) * clip_length + clip_length,LengthFiles); 
    unary_res = cell(length(frame_range),1);
    binary_res = cell(length(frame_range) - 1,1);
    att_ratio_res = cell(length(frame_range),1);
    att_level_res = cell(length(frame_range),1);
        
    for frame_clip_i = 1 : length(frame_range)
        
        frame = frame_range(frame_clip_i);
        subjects_gaze = all_sub_gaze{frame};
        
%         img_path = [filePath,'/Images/','V0-S1-G_1','/frame_sel/'];
%         params.framenum = num2str(frame,'%04d');
%         img_top_name = [params.framenum,'.jpg'];
%         data.img_top = imread([img_path,'top/',img_top_name]);

        focus_level = cell(size(subjects_gaze,1),3);
        for can_i = 1:length(subjects_gaze)
            candidate =  subjects_gaze{can_i};
            focus_level_i = single_focus(params,candidate,subjects_gaze,can_i);
            focus_level(can_i,:) = focus_level_i;
        end

        attention = cell(size(focus_level,1),1);
        att_level = zeros(size(focus_level,1),1);
        att_ratio = zeros(size(focus_level,1),1);
        
        for i = 1 : size(focus_level,1)
            att_level(i) = cell2mat(focus_level(i,1));
            att.level = att_level(i);
            att_ratio(i) = cell2mat(focus_level(i,2))/size(subjects_gaze,1);
            att.ratio =  att_ratio(i);
            % att.score = -log(0.5 * att.level + 0.5 * att.ratio);
            attention{i} = att;
        end

        %  the spatialtemporal smooth of frame t and frame t+1 
        if frame_clip_i > 1

            last_subjects_gaze = data.last_subjects_gaze;
            last_attention = data.last_attention;

            d_att_level = zeros(size(last_subjects_gaze,1),size(subjects_gaze,1));
            d_att_ratio = zeros(size(last_subjects_gaze,1),size(subjects_gaze,1));
            overlap_t_r = zeros(size(last_subjects_gaze,1),size(subjects_gaze,1));

            for can_t_i = 1:size(last_subjects_gaze,1)  %  All candidates from frame_t
                candidate_t_i =  last_subjects_gaze{can_t_i};
                attention_t_i = last_attention{can_t_i};

                for can_r_i = 1:size(subjects_gaze,1)   %  All candidates from frame_r
                    candidate_r_i =  subjects_gaze{can_r_i};
                    attention_r_i = attention{can_r_i};

                    d_att_level(can_t_i,can_r_i) = attention_t_i.level - attention_r_i.level;
                    d_att_ratio(can_t_i,can_r_i) = attention_t_i.ratio - attention_r_i.ratio;
                    overlap_t_r(can_t_i,can_r_i) = get_IOU(candidate_t_i.objs_top,candidate_r_i.objs_top);

                end
            end
        end

        data.last_subjects_gaze = subjects_gaze;
        data.last_attention = attention;
        
        weight_att_ratio = 1;
        weight_att_level = 1;
        weight_att_change = 1;
        weight_att_overlap = 1;

        %unary_res{frame_clip_i} =  -(att_level + att_ratio);
        unary_res{frame_clip_i} = -1./(1+exp(-(weight_att_level * att_level + weight_att_ratio * att_ratio)));
        att_ratio_res{frame_clip_i} = att_ratio;
        att_level_res{frame_clip_i} = att_level;
        
        if frame_clip_i > 1
            %binary_res{frame_clip_i-1} = -(overlap_t_r + d_att_level + d_att_ratio);
            binary_res{frame_clip_i-1} = -1./(1+exp(-(weight_att_overlap * overlap_t_r + weight_att_change * (d_att_level + d_att_ratio))));
        end

    end 
    unary_res;
    binary_res;
    att_ratio_res;
    T_frame = length(frame_range);
    N = 10;

    for i = 1:(T_frame-1)
        for j = 1:(length(unary_res{i}))
            A{i}(j,:) = binary_res{i}(j,:) + unary_res{i}(j); %%
        end
    end
    
    A{T_frame} = unary_res{T_frame,1};

    candidates = zeros(T_frame-1,1);
    maxMat = zeros(N,T_frame);

    maxRow = zeros(N,1);

    sum_scores = 0;
    
    for t = T_frame:-1:1
        maxRow = maxMat(:,1);
        for i = 1:size(A{t})
            [a,b] = min(A{t}(i,:)'+maxRow(1:size(A{t},2)));
    %         maxMat
    %         A{t,1}(i,:)'
    %         maxRow(1:size(A{t,1},2))+        A{t,1}(i,:)'
            maxMat(i,1) = a;
            maxMat(i,t+1) = b;
    %         [maxMat(i,T+1),maxMat(i,t-1)] = [a,b];
    %         [a,b] = min(A{t,1}(i,:)'+max_mat(1:size(A{t,1},2),T+1))
    %         (A{t,1}(i,:))'+ max_mat(1:size(A{t,1},2),T+1)
    %          max_mat(1:size(A{t,1},2),T+1)
    %          A{t,1}(i,:)
        end
    end
    
    for t = 1:T_frame
        if t == 1
    %     maxScore = maxMat(1,1);
            [~,candidates(t)] = min(unary_res{t});
    %     elseif t == T
    %         [~,candidates(t)] = min(AA);
        else
            candidates(t) = maxMat(candidates(t-1),t);
        end
    end
    
    all_att_ratio = zeros(length(att_ratio_res),1);
    all_att_level = zeros(length(att_ratio_res),1);
    for i = 1 : length(att_ratio_res)
        att_ratio_i = att_ratio_res{i};
        att_level_i = att_level_res{i};
        candidates_i = candidates(i);
        all_att_ratio(i) = att_ratio_i(candidates_i);
        all_att_level(i) = att_level_i(candidates_i);
    end
    avg_att_ratio = mean(all_att_ratio);
    avg_att_level = mean(all_att_level); 
    
%   avg_all = 1./(1+exp(-(avg_att_ratio + avg_att_level)));
    avg_all = avg_att_ratio + avg_att_level;
%   if avg_att_ratio < 0.45

    threshold = 1;
    if params.motion_dir == 1
        threshold = 0.3;   % decrease the threshold when using motion direction
    end
    if avg_all < threshold
        candidates(:) = 0;
    end
    
%     for t = 1:T_frame
%         if t == 1
%             meanCandidates = mean(unary_res{t},2);
%             outTest = filloutliers(meanCandidates,0,'grubbs');
%             if (norm(outTest - meanCandidates)) == 0
%                 candidates(t) = 0;
%             end
%         else
% %             meanCandidates = mean(A{t},2);
% %             outTestPro = filloutliers(meanCandidates,0,'median');
% % %             outTestPro = filloutliers(meanCandidates,0,'percentiles',[10,100]);
% %             if (norm(outTestPro - meanCandidates)) == 0
% %                 candidates(t) = 0;
% %             end
%           %AAAA = [mean(mean(A{t},2)),median(mean(A{t},2)),var(mean(A{t},2)),(clip_i-1)*50+t]
% %         if min(mean(A{t}))<-1.3
% %             candidates(t) = 0;
% %         end
% % 
% %         if (norm(outTest-outTestPro)) ==0
% %             candidates(t) = 0;
% %         end
% %             meanCandidates = mean(binary_res{t-1},2);
% %             outTestPro = filloutliers(meanCandidates,0,'grubbs','ThresholdFactor',0.95);
% %             %outTestPro = filloutliers(meanCandidates,0,'percentiles',[7,100]);
% %             if (norm(outTestPro - meanCandidates)) == 0
% %                 candidates(t) = 0;
% %             end
%         end
%     end

        %counter0 = 0;
        setC = zeros(2,length(frame_range));
        for j = 1:length(frame_range)
            jj = clip_length*(clip_i-1)+j;
            disp(['Frame:' num2str(jj)]);
            gaze_frm = all_sub_gaze{jj};
%             if candidates(j) == 0 
%                 counter0 = counter0 + 1;
%             else 
            if candidates(j) ~= 0 && candidates(j)<=length(gaze_frm)
                gaze_sub = gaze_frm{candidates(j)};
                center_x = gaze_sub.objs_top(1) + 0.5 * gaze_sub.objs_top(3);
                center_y = gaze_sub.objs_top(2) + 0.5 * gaze_sub.objs_top(4);
                setC(:,j) = [center_x;center_y];                
            end
        end
        clusterK = 3;
        clusterC = clusterdata(setC',clusterK);
        mostCandidates = mode(clusterC);
%         if counter0 > mostCandidates
%             bestC(i) = 0;
%         else
            indexC = find(clusterC==mostCandidates);

            bestC(clip_i) = candidates(indexC(1));

            if bestC(clip_i)>0

                if bestC(clip_i) <= length(gaze_frm)
                gaze_sub = gaze_frm{bestC(clip_i)};
                end

                center_x = gaze_sub.objs_top(1) + 0.5 * gaze_sub.objs_top(3);
                center_y = gaze_sub.objs_top(2) + 0.5 * gaze_sub.objs_top(4);
                candidateOfClip(:,clip_i) = [center_x;center_y];
            else
                candidateOfClip(:,clip_i) = [0;0];
            end
%         end
        clip_res = candidates;  % Put your CAP detection results here
        video_res(frame_range) = clip_res';
 
        clear A;
    end

    %% Save the CAP detection results   
    save(fullfile(params.data_directory,'Data',['CAP_' params.scene_i.scene '.mat']),'video_res','-v7.3');
    if params.motion_dir == 1
        save(fullfile(params.data_directory,'Motion',['CAP_' params.scene_i.scene '.mat']),'video_res','-v7.3');
    end
else
    
    load(fullfile(params.data_directory,'Data',['CAP_' params.scene_i.scene '.mat']));

end
   
%% Visualize the gaze and CAP results
   
if params.vis_CAP_res == 1
    for frame = 1 : size(all_sub_gaze,1)
        
        figure(1)
        framenum = num2str(frame,'%04d');
        img_top_name = [framenum,'.jpg'];
        img_top = imread([img_path(1:end-2),'T','/top/',img_top_name]);
        imshow(img_top);hold on; 
        
        gaze_frm = all_sub_gaze{frame};
        text(100,100,num2str(framenum),'color','y','fontsize',25);
        
        for sub_i =  1 : size(gaze_frm,1)
            
            gaze_sub = gaze_frm{sub_i};
            
            center_x = gaze_sub.objs_top(1) + 0.5 * gaze_sub.objs_top(3);
            center_y = gaze_sub.objs_top(2) + 0.5 * gaze_sub.objs_top(4);
            gaze_line_h = 50;
            gaze_x = center_x + gaze_line_h * gaze_sub.view_angle(1);
            gaze_y = center_y - gaze_line_h * gaze_sub.view_angle(2);

            plot([center_x, gaze_x], [center_y, gaze_y],'Linewidth',1.5,'color',[1,0.5,0]);
            text(center_x-30,center_y-30,num2str(sub_i),'color','g','fontsize',15);
 
        end
        
        CAP_id = video_res(frame);
         
        if  CAP_id > 0
        CAP = gaze_frm{CAP_id};
            
        center_x = CAP.objs_top(1) + 0.5 * CAP.objs_top(3);
        center_y = CAP.objs_top(2) + 0.5 * CAP.objs_top(4);

        plot(center_x,center_y,'rp','MarkerSize',15,'MarkerFaceColor','r');
        else
            
        plot(175,150,'rp','MarkerSize',15,'MarkerFaceColor','r');
        end
        
%       clip_i = floor(frame/clip_length)+1;
%       center_x = candidateOfClip(1,clip_i);
%       center_y = candidateOfClip(2,clip_i);
%       plot(center_x,center_y,'rp','MarkerSize',15,'MarkerFaceColor','r');        
        hold off;
%       close off;
        pause(0.05);
  
    end
    
end
  
end





