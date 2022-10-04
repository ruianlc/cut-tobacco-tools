clc;clear;close all
%%
% --- 说明 ---
% 读取多个子目录下的片烟图像
% 对每张图像中的片烟进行模拟切丝，获取图像中每个片烟分割出来的烟丝长度
% 
% Programmer: Robin An, 2021.07.28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO
% DONE
% 1.多目录读取
% 2.count_len函数优化
% 3.切割产生断丝
% 4.不填充，片烟内部存在空洞。但一旦考虑此点，需要判断每个片烟内部是否存在闭合区域，即片烟所在的闭合区域填充，但片烟内部存在的闭合区域不做填充
% 5.剔除异常原始数据
% ===
% 6.最小外接圆、最大内切圆
% 7.最小外接矩

rootdir = 'G:\workspace\数据存储\制丝段离线实验\片烟检测-2020.07.21';
% rootdir = 'G:\workspace\数据存储\制丝段离线实验\少量图片测试';
subdir = dir(rootdir);

% 大中小碎分别设定不同阈值
big_threshold = 800;
middle_threshold = 500;
small_threshold = 200;
detritus_threshold = 100;

% 面积
% K = 0.0175;
% B = 0;

% 面积
K = 0.0241;
B = 0;

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
    Atotal = cell(num_img,1); % cell存储每张图片烟丝的面积
    Ntotal = cell(num_img,1); % cell存储每张图片的名称
    
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
        
       %% 二值图像    
        gray = rgb2gray(I); % 图像的灰度处理
        level = graythresh(gray);
        bw = imbinarize(gray,level);  % 图像的二值化处理
        bw1 = ~bw; % 获得每个片烟
        bw2 =imfill(bw1,'holes');
        adjustImage = bwareaopen(bw2, threshold); % 剔除小面积图像
       
       %% 提取每个连通区域和相应的最小外接圆
        regions = regionprops(adjustImage); % 计算联通区域，region = [{area,centroid,boundingbox}]
        [boundaries,L] = bwboundaries(adjustImage); % 获取边界，此处需对bw1做空洞填充，使边界数量与连通区域一致
        
        [label,sheet_num] = bwlabel(adjustImage);
        Images = cell(1,sheet_num); % 每个片烟所在区域
        tarea = cell(1,sheet_num); % 每个片烟面积
        
        tcenter = cell(1,length(boundaries));
        tradius = zeros(1,length(boundaries));
        circlearea = zeros(1,length(boundaries));
        rrate = zeros(1,length(boundaries));
        
        tr = cell(1,sheet_num);
        for k = 1:sheet_num % 对每个片烟分别操作
            tarea{k} = regions(k).Area;
            if tarea{k} > 1000000
                tarea{k} = []; %最大片烟84万，设定片烟阈值不超过100万
                continue;
            end
            boundary = boundaries{k};
            [center,radius] = minboundcircle(boundary(:,2),boundary(:,1)); % 最小外接圆半径和圆心
            tcenter{k} = center;
            tradius(k) = radius;
            circlearea(k) = pi * (radius.^2);
            
            rrate(k) = tarea{k} / circlearea(k);
            tr{k} = [tarea{k},rrate(k)];
        end
        Ntotal{j} = subdirfiles(j).name(1:end-4); 
        Atotal{j} = [tr{:}];
    end % j - 每个目录下的图像
    
    %% 统计结果：片烟面积、外接圆面积、圆度率
    [curpath,~] = fileparts(mfilename('fullpath')); % 当前目录
    dsave = strcat(curpath,'\片烟面积统计.xlsx');
    writecell(Ntotal,dsave,'Sheet',i-2,'Range','A1'); % 第一列写入图像文件名
    writecell(Atotal,dsave,'Sheet',i-2,'Range','B1'); % 写入片烟面积
end %i - 子目录