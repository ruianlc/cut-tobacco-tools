clc;clear;close all
%%
% --- 说明 ---
% 对片烟图像进行模拟切丝，获取图像中每个片烟分割出来的烟丝长度
% 
% Programmer: Robin An, 2021.07.27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TODO 
% 最小外接矩、最大内切圆、最小外切圆、连通区域切割、不填充

 pathfile = 'G:\workspace\数据存储\制丝段离线实验\片烟检测-2020.07.21\大片\20210721151138.bmp';
% pathfile = 'G:\workspace\数据存储\制丝段离线实验\片烟检测-2020.07.21\中片\20210721162841.bmp';
% pathfile = 'G:\workspace\数据存储\制丝段离线实验\20210721162841.bmp';
% pathfile = 'G:\workspace\数据存储\制丝段离线实验\片烟检测-2020.07.21\小片\20210721195251.bmp';
% pathfile = 'G:\workspace\数据存储\制丝段离线实验\片烟检测-2020.07.21\碎片\20210721230318.bmp';

I = imread(pathfile);
figure
imshow(I)
title('原始图像')
%I = source(120:end,:,:);
%imtool(I)

%% 二值图像
rgb = reshape(I,size(I,1)*size(I,2),3);
rgb(rgb(:,3)>150,:) = 0; % 去除背景
%rgb(rgb(:,1)<30,:) = 0;  % 去除影响元素

rgb(rgb(:,3)~= 0,:) = 1;
img = reshape(rgb,size(I,1),size(I,2),3);
bw = img(:,:,1);
BW = logical(bw);
figure
imshow(BW)
title('二值化图像')

% se = strel('disk',5); % 这里是创建一个半径为5的平坦型圆盘结构元素
% BW2 = imerode(BW1,se);% 腐蚀
% figure
% imshow(BW2)
% title('二值化腐蚀图像')
% 
% montage({I,BW,BW1,BW2},'Size',[2 2])
adjustImage = bwareaopen(BW, 300); % 剔除小面积图像
figure
imshow(adjustImage)
title('图像二值化剔除小面积')

%% 获取区域边界
[boundaries,L] = bwboundaries(adjustImage);
% tcenter = cell(1,length(boundaries));
% tradius = zeros(1,length(boundaries));
area = zeros(1,length(boundaries));  % 所有片烟面积
tarea = zeros(1,length(boundaries)); % 所有片烟外接圆面积
rrate = zeros(1,length(boundaries)); % 所有片烟圆度率

imshow(label2rgb(L, @jet, [.5 .5 .5]))
hold on
for k = 1:length(boundaries)
    boundary = boundaries{k};
    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end

%% 提取每个连通区域
regions = regionprops(adjustImage); % 计算连通区域，region = [{area,centroid,boundingbox}]
[label,sheet_num] = bwlabel(adjustImage);

% 方框标注连通区域
tpos=regionprops(label,'BoundingBox');
tcentroid = regionprops(label,'Centroid');
for i = 1:sheet_num % 每个片烟
    rectangle('position',tpos(i).BoundingBox,'edgecolor','r', 'LineWidth',1.5);
    line([tcentroid(i,1).Centroid(1,1), tcentroid(i,1).Centroid(1,1)],[tcentroid(i,1).Centroid(1,2)-35, tcentroid(i,1).Centroid(1,2)+35],'Color','b', 'LineWidth',2);
    text(tcentroid(i,1).Centroid(1,1),tcentroid(i,1).Centroid(1,2), num2str(i),'Color', 'r');
end

K = 0.1314;
B = 0;
% K = 0.1336;
% B = -0.6188;

PL = 13; % 像素与实际长度，切丝宽度1mm，对应像素为12.1168；切丝宽度2mm，对应像素为19
Images = cell(1,sheet_num); % 每个图像所有片烟所在区域

LL = cell(sheet_num,1);% 每张图像片烟指标：面积、圆度率、烟丝长度

[cur_path,~] = fileparts(mfilename('fullpath'));

for i = 1:sheet_num % 对每个片烟分别操作
    bw = label;
    pos = regions(i).BoundingBox;
    
    r1 = round(pos(2));
    c1 = round(pos(1));
    w = pos(3);
    h = pos(4);
    r2 = r1+h-1;
    c2 = c1+w-1;
    bw(bw ~= i) = 0;  % 其他像素置0
    bw(bw == i) = 1;  % 将等于i的像素置1
    Images{i} = bw(r1:r2,c1:c2);  % 每个片烟的像素区域
    nImage = ~Images{i};
        
    %% 片烟切割
    stick_num = floor(w/PL);
    
    %% 片烟切割烟丝长度
    savepic = strcat(cur_path,'\Figure\',num2str(i),'.png');
    h0 = figure('color','w');
    set(h0,'Visible','off');
    
    for j = 1:stick_num-1 
        nImage(:,j*PL+1:j*PL+2) = 1;
    end
    imshow(nImage);
    print(h0,'-dpng',savepic);
end
