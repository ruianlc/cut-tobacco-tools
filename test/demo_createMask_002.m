clc;clear;close all
fprintf('starts...')

ipath = 'C:\Users\Administrator\Documents\MATLAB\cut-tobacco-tools\data\002.bmp';
imgData = imread(ipath);

% 切割边缘长度
cuttingLenUp = 200;    % 上边缘
cuttingLenBot = 200;  % 下边缘 
cuttingLenLeft = 200;  % 左边缘
cuttingLenRight = 200;% 右边缘

% 去除小面积阈值
areaThreshold = 500;

[BW,~] = createMask_002(imgData);
BW1 = BW(cuttingLenUp:end-cuttingLenBot,cuttingLenLeft:end-cuttingLenRight);  % 图像切割
BW2 = ~BW1; % 二值化图像
figure
imshow(BW2)

adjustImage = bwareaopen(BW2, areaThreshold);    % 剔除小面积图像
figure
imshow(adjustImage)

regions = regionprops(adjustImage); % 计算连通区域，region = [{area,centroid,boundingbox}]
[label,sheet_num] = bwlabel(adjustImage);

% 方框标注连通区域
tpos=regionprops(label,'BoundingBox');
tcentroid = regionprops(label,'Centroid');
for i = 1:sheet_num
    rectangle('position',tpos(i).BoundingBox,'edgecolor','r', 'LineWidth',1.5);
    line([tcentroid(i,1).Centroid(1,1), tcentroid(i,1).Centroid(1,1)],[tcentroid(i,1).Centroid(1,2)-35, tcentroid(i,1).Centroid(1,2)+35],'Color','b', 'LineWidth',2);
    text(tcentroid(i,1).Centroid(1,1),tcentroid(i,1).Centroid(1,2), num2str(i),'Color', 'r');
end


%% 计算联通区域

Lens = zeros(sheet_num, 1);
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
    Lens(i) = sum(DDist);
end


fprintf('the end...')