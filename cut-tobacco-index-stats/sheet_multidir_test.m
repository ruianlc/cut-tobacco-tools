clc;clear;close all
%%
% --- 说明 ---
% 读取多个子目录下的片烟图像
% 对每张图像中的片烟进行模拟切丝，获取图像中每个片烟分割出来的烟丝长度
% 
% Programmer: Robin An, 2021.07.30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO( DONE now )
% 1.多目录读取
% 2.count_len函数优化
% 3.切割产生断丝
% 4.不填充，片烟内部存在空洞。但一旦考虑此点，需要判断每个片烟内部是否存在闭合区域，即片烟所在的闭合区域填充，但片烟内部存在的闭合区域不做填充
% 5.剔除异常原始数据
% 6.最小外接圆
% 7.最小外接矩

%rootdir = 'D:\workspace\datasource\上烟制丝卷包\在线制丝实验\A';
% rootdir = 'D:\workspace\datasource\上烟制丝卷包\在线制丝实验\B';
 rootdir = 'D:\workspace\datasource\上烟制丝卷包\在线制丝实验\B小样';

subdir = dir(rootdir);

% 大中小碎分别设定不同阈值
big_threshold = 800;
middle_threshold = 500;
small_threshold = 200;
detritus_threshold = 100;

for i = 1:length(subdir)
    
    % 跳过'.'和'..'目录
    if(strcmp(subdir(i).name,'.')||...
            strcmp(subdir(i).name,'..')||...
            ~subdir(i).isdir())
        continue;
    end
    
    % 子文件夹下找后缀为bmp的文件
    subdirpath = fullfile( rootdir, subdir(i).name);
    subdirfiles = dir(strcat(subdirpath,'\','*.bmp'));
    num_img = length(subdirfiles); % 每个子文件夹图像总数
    Ltotal = {};%cell(num_img,1); % cell存储每张图片烟丝的面积
    
    if(isequal(num_img,0))
        return;
    end
    
    if contains(subdir(i).name,'大')
        threshold = big_threshold;
    elseif contains(subdir(i).name,'中')
        threshold = middle_threshold;
    elseif contains(subdir(i).name,'小')
        threshold = small_threshold;
    elseif contains(subdir(i).name,'碎')
        threshold = detritus_threshold;
    else
        threshold = middle_threshold;
    end
        
    for j = 1:num_img
        
        imgfile_path = fullfile(subdirpath, subdirfiles(j).name);
        I = imread(imgfile_path);
        
%         rgb = reshape(I,size(I,1)*size(I,2),3);
%         rgb(rgb(:,3)>110,:) = 0; % 根据B值差异去除背景
%         rgb(rgb(:,3)~= 0,:) = 1;
%         img = reshape(rgb,size(I,1),size(I,2),3);
%         bw = img(:,:,1);
%         BW = logical(bw);  
%         adjustImage = bwareaopen(BW, threshold); % 剔除小面积图像
        
        gray = rgb2gray(I); % 图像的灰度处理
        level = graythresh(gray);
        bw = imbinarize(gray,level);  % 图像的二值化处理
        bw1 = ~bw; % 获得每个片烟
        bw2 = imfill(bw1,'holes');
        adjustImage = bwareaopen(bw2, threshold); % 剔除小面积图像

        figure
        imshow(adjustImage)
    end % j - 每个目录下的图像
    
end %i - 子目录