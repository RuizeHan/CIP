function [objs_top,objs_ego] = frame2obj(filePath,params)
   
%     tv_mat_name = [params.scenenum,'_top.mat'];
%     hv_mat_name = [params.scenenum,'_hor.mat'];
%     s = whos('-file',[filePath,'/GT/',tv_mat_name]);
%     tv_GT = load([filePath,'/GT/',tv_mat_name]);
%     hv_GT = load([filePath,'/GT/',hv_mat_name]);   
%     
%     anno_top = cell2mat(tv_GT.(s.name)(params.frame_i,2));
%     anno_ego = cell2mat(hv_GT.(s.name)(params.frame_i,2));
%     
%     anno_top(:,3) = anno_top(:,3) - anno_top(:,1);
%     anno_top(:,4) = anno_top(:,4) - anno_top(:,2);
%     anno_ego(:,3) = anno_ego(:,3) - anno_ego(:,1);
%     anno_ego(:,4) = anno_ego(:,4) - anno_ego(:,2);
    
    % translate the GT into MHA formate: -1:Gopro,0:unmatched object
%     tmp_lab_top = anno_top(:,5);
%     tmp_lab_ego = anno_ego(:,5);
%     gopro_idx = find(tmp_lab_top==0);
%     tv_isnum_idx = ismember(tmp_lab_top,tmp_lab_ego);
%     hv_isnum_idx = ismember(tmp_lab_ego,tmp_lab_top);
%     tmp_lab_ego(hv_isnum_idx==0) = 0;
%     tmp_lab_top(tv_isnum_idx==0) = 0;
%     tmp_lab_top(gopro_idx) = -1;
%     anno_top(:,5) = tmp_lab_top;
%     anno_ego(:,5) = tmp_lab_ego;
    
    % load the data as MHA formate
%     cam_top = anno_top((anno_top(:,5)== -1),1:4);   
%     lab.top = anno_top((anno_top(:,5)~= -1),5);    
%     lab.ego = anno_ego(:,5);
%     objs_ego_gt = anno_ego(:,1:4);
%     objs_top_gt = anno_top((anno_top(:,5)~= -1),1:4);
    
%     cam_search = anno_top(:,1:4);  %all the objs in the top view to search for the ego-camera
    
    if params.detection_method ~= 1
        det_filePath = [filePath,'Detection/'];
        Det = load([det_filePath,params.scenenum,'_hor.mat']);
        Det_tv = load([det_filePath,params.scenenum,'_top.mat']);
        s = whos('-file',[det_filePath,params.scenenum,'_hor.mat']);
        det_res = cell2mat(Det.(s.name)(params.frame,2));
        objs_can = [det_res(:,3)-det_res(:,1),det_res(:,4)-det_res(:,2)];       
        det_d_res = cell2mat(Det_tv.(s.name)(params.frame,2));

        if size(det_d_res,1) ~=0
            det_top(:,1:2) = det_d_res(:,1:2);
            det_top(:,3:4) = [det_d_res(:,3)-det_d_res(:,1),det_d_res(:,4)-det_d_res(:,2)];   
        else
            det_top = [1,1,1,1];
        end
%         cen_dis = zeros(size(det_top,1));
%         for i = 1:size(det_top,1)
%             cen_dis(i) = get_center_dis(cam_top,det_top(i,:));
%         end            
    end
    
switch  params.detection_method
    
    case 1 % GT detection           
        objs_ego = objs_ego_gt;
        objs_top = objs_top_gt;
        if isempty(objs_top_gt)
            objs_top = [1,1,0,0];
        end
                
    case 2   % raw detection results        
        obj_idx = find(objs_can(:,2) > 0);   %remove the bounding box smaller than the minimum
        objs_ego(:,1:2) = det_res(obj_idx,1:2);
        objs_ego(:,3) = objs_can(obj_idx,1);
        objs_ego(:,4) = objs_can(obj_idx,2);
        objs_ego(:,5) = det_res(obj_idx,5);
        objs_ego(:,6) = det_res(obj_idx,6);
        % objs_top = det_top((cen_dis>20),:);   % get the top detection results except for the ego-camera wearer
        objs_top = det_top;   % get all the top detection results
                        
    case 3 % detection with GT selection 
               
        objs_top = det_top((cen_dis>20),:);   % get the top detection results except for the ego-camera wearer
        objs_ego(:,1:2) = det_res(:,1:2);
        objs_ego(:,3) = objs_can(:,1);
        objs_ego(:,4) = objs_can(:,2);
                
        top_dis = zeros(size(objs_top_gt,1),size(objs_top,1));
        ego_IOU = zeros(size(objs_ego_gt,1),size(objs_ego,1));
        dis_thr_top = 20;
        for i = 1 : size(objs_top_gt,1)
            for j = 1 : size(objs_top,1)  
                top_dis(i,j) = get_center_dis(objs_top_gt(i,:),objs_top(j,:)); 
            end
        end
        [min_col,~] = min(top_dis,[],1); 
        objs_top = objs_top((min_col < dis_thr_top),:);
        if isempty(objs_top)
            objs_top = [1,1,0,0];
        end
        
        IOU_thr = 0.5;
        for m = 1 : size(objs_ego_gt,1)
            for n = 1 : size(objs_ego,1)            
             ego_IOU(m,n) = get_IOU(objs_ego_gt(m,:),objs_ego(n,:));   
            end
        end
        [max_col,~] = max(ego_IOU,[],1);
        objs_ego = objs_ego((max_col > IOU_thr),:);              
end

end