fileName = 'D:\Research\multi_camera\aligned\Multi-video_correlation\scene1_55\g_3_0.mp4';
str = 'D:\Research\multi_camera\video2jpg\GoPro3\';
obj = VideoReader(fileName);
numFrames = obj.NumberOfFrames;% 帧的总数

for k = 1 : numFrames% 读取数据
    frame = read(obj,k);
    %imshow(frame);%显示帧
    c='0000';
    b= num2str(k);
    b=[c(1:4-length(b)) b];
    imwrite(frame,strcat(str,b,'.jpg'),'jpg');% 保存帧
end
