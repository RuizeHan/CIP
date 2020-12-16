function params = params_setting(params)

params.search_gopro = 2; %  0:GT 1: auto,  2:Tracking res
params.detection_method = 2; % 1: GT, 2:Raw Detection  3: Detection with selection
params.Rho = 25; % plenty factor
params.x_thresh = 100; 
params.dis_thresh = 1000; 
params.Lambda = 0.015; % weight for y defined euclidean distance 0.05
params.occ = 1;
params.occ_angle = 2;
params.fast_search = 0;
params.search_ite = 6;
params.interval = 1;
% params.search_num = ceil(360/params.interval);
params.temp_angle = 0;
params.min_scl_degree = 0.1;   % the min degree search range
params.gma = 0.05; % Motion explotion Average parameter
params.num_gopro = 3;

end