fileName = 'D:\Research\multi_camera\aligned\Multi-video_correlation\scene1_55\g_3_0.mp4';
str = 'D:\Research\multi_camera\video2jpg\GoPro3\';
obj = VideoReader(fileName);
numFrames = obj.NumberOfFrames;% ֡������

for k = 1 : numFrames% ��ȡ����
    frame = read(obj,k);
    %imshow(frame);%��ʾ֡
    c='0000';
    b= num2str(k);
    b=[c(1:4-length(b)) b];
    imwrite(frame,strcat(str,b,'.jpg'),'jpg');% ����֡
end
