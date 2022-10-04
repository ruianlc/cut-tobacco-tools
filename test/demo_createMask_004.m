clc;clear;close all
fprintf('starts...')

ipath = 'C:\Users\Administrator\Documents\MATLAB\cut-tobacco-tools\data\004.bmp';
imgData = imread(ipath);

% 切割边缘长度
cuttingLenUp = 20;    % 上边缘
cuttingLenBot = 20;  % 下边缘 
cuttingLenLeft = 20;  % 左边缘
cuttingLenRight = 20;% 右边缘

% 去除小面积阈值
areaThreshold = 500;
resiThreshold = 1;

[BW,~] = createMask_004(imgData);
BW1 = BW(cuttingLenUp:end-cuttingLenBot,cuttingLenLeft:end-cuttingLenRight);  % 图像切割
BW2 = ~BW1; % 二值化图像

adjustImage = bwareaopen(BW2, areaThreshold);    % 剔除小面积图像
% figure
% imshow(adjustImage)

%% 计算联通区域
regions = regionprops(adjustImage); % 计算连通区域，regions = [{area,centroid,boundingbox}]
[label,sheet_num] = bwlabel(adjustImage);

Lens = zeros(sheet_num, 1);
lineFlags = zeros(sheet_num,1);
residualsSums = zeros(sheet_num,1);

for i = 1:sheet_num
    bw = label;
    pos = regions(i).BoundingBox;
    r1 = round(pos(2));
    c1 = round(pos(1));
    w = pos(3);
    h = pos(4);
    r2 = r1+h-1;
    c2 = c1+w-1;
    bw(bw ~= i) = 0;  % 其他像素置
    bw(bw == i) = 1;  % 将等于i的像素置1
	subImg = bw(r1:r2,c1:c2,:);  % 单个烟丝二值化图像
  
    skL = sketelon(subImg); % 提取骨架
    [rows,columns] = find(skL == 1);  % 提取骨架点坐标
    DDist = count_len(rows,columns);
    [lineFlag,residualsSum] = check_on_line(rows,columns,resiThreshold);
    
    lineFlags(i) = lineFlag;
    residualsSums(i) = residualsSum;
    Lens(i) = sum(DDist);
end


fprintf('the end...')