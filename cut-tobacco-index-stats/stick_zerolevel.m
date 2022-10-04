clc;clear;close all
%%
% [filename,pathname,~] = uigetfile({'*.bmp';'*.jpg'},'选择图像');  % 选择文件夹
% pathfile = strcat(pathname,'\',filename);  % 组合路径+文件名
% I = imread(char(pathfile));  % 读取图片
%pathfile = 'G:\workspace\数据存储\上烟制丝卷包\在线实验-烟丝长度检测\#01\1.长\2021-09-07_11_01_00_041.bmp';
%pathfile = 'G:\workspace\数据存储\上烟制丝卷包\在线实验-烟丝长度检测\#01\1.长\2021-09-07_11_22_25_344.bmp';
%pathfile = 'G:\workspace\数据存储\上烟制丝卷包\在线实验-烟丝长度检测\#02\1.长\2021-09-07_13_06_04_216.bmp';
%pathfile = 'G:\workspace\数据存储\上烟制丝卷包\在线实验-烟丝长度检测\#01\2.中\2021-09-07_11_45_26_790.bmp';
%pathfile = 'G:\workspace\数据存储\上烟制丝卷包\在线实验-烟丝长度检测\#04\2.中\2021-09-07_14_26_28_168.bmp';
%pathfile = 'G:\workspace\数据存储\上烟制丝卷包\在线实验-烟丝长度检测\#06\3.短\2021-09-07_15_34_46_247.bmp';
%pathfile = 'G:\workspace\数据存储\上烟制丝卷包\在线实验-烟丝长度检测\#01\3.短\2021-09-07_11_43_36_383.bmp';
%pathfile = 'G:\workspace\数据存储\上烟制丝卷包\在线实验-烟丝长度检测\#01\4.碎\2021-09-07_12_00_02_420.bmp';
%pathfile = 'G:\workspace\数据存储\上烟制丝卷包\在线实验-烟丝长度检测\#01\4.碎\2021-09-07_13_03_35_649.bmp';

%pathfile = 'G:\workspace\数据存储\上烟制丝卷包\在线实验-烟丝长度检测\#06\1.长\2021-09-07_15_39_52_775.bmp';
%pathfile = 'D:\workspace\科技项目\上海烟草\卷接-制丝项目\上烟制丝卷包\实验方案\卷接段实验\烟支拆开摆放烟丝_1支12-15分钟 - 副本.jpg';
pathfile = 'D:\workspace\科技项目\上海烟草\卷接-制丝项目\上烟制丝卷包\实验方案\卷接段实验\烟支拆开摆放烟丝_1支12-15分钟 - 副本 (2).jpg';
threshold = 1;

I = imread(pathfile);
figure
imshow(I)
%             titleString = strcat('第',num2str(i),'个样本 -- ','文件夹：',subsubdir(j).name,' -- 第',num2str(s),'张图像');
%             title(titleString)

%             K = 0.3262;
%             B = 0.7791;
%             K = 0.1314; % 标样尺寸
%             B = 0;
%             pixel = 120;
%             L = K * pixel + B
%
%             100 13.1
%             120 15.8
%             150 19.7
%             180 23.7
%             200 26.3

source = I;
% source = I(1:end-100,200:end-200,:); %截去左下角阴影影响
% figure
% imshow(source)
% title('处理图像')

rgb = reshape(source,size(source,1)*size(source,2),3);
rgb(rgb(:,3) > 200,:) = 0;
rgb(rgb(:,3)~= 0,:) = 1;
img = reshape(rgb,size(source,1),size(source,2),3);
bw0 = img(:,:,1);
bw = logical(bw0);


figure
imshow(bw)
title('二值化图像')

BW2 = imfill(bw,'holes');  % 填充图像
figure
imshow(BW2)
title('二值化填充图像')
regions0 = regionprops(BW2);

adjustImage = bwareaopen(BW2, threshold); % 剔除小面积图像
figure
imshow(adjustImage)
title('二值化剔除小面积图像')

regions = regionprops(adjustImage); % 计算联通区域
[labels,Num] = bwlabel(adjustImage);

% tpos=regionprops(labels,'BoundingBox');
% tcentroid = regionprops(labels,'Centroid');
% for i = 1:Num % 每个片烟
%     rectangle('position',tpos(i).BoundingBox,'edgecolor','r', 'LineWidth',1.5);
%     line([tcentroid(i,1).Centroid(1,1), tcentroid(i,1).Centroid(1,1)],[tcentroid(i,1).Centroid(1,2)-35, tcentroid(i,1).Centroid(1,2)+35],'Color','b', 'LineWidth',2);
%     text(tcentroid(i,1).Centroid(1,1),tcentroid(i,1).Centroid(1,2), num2str(i),'Color', 'r');
% end

% 提取根个烟丝图像
% K = 0.3262;
% B = 0.7791;
K = 1;
B = 0;
Lens = zeros(1,Num); % 单张图像中所有烟丝的长度

for t = 1:Num % idx为单张图像中每个烟丝的序号
    tlabel = labels;
    
    pos = regions(t).BoundingBox; % 最小外接矩形的x，y边界（下界和区间），pos = [x,y,width,height]
    r1 = round(pos(2));
    c1 = round(pos(1));
    w = pos(3);
    h = pos(4);
    r2 = r1+h-1;
    c2 = c1+w-1;
    tlabel(tlabel ~= t) = 0;  % 其他像素置
    tlabel(tlabel == t) = 1;  % 将等于i的像素置1
    stick = tlabel(r1:r2,c1:c2);  % 单个烟丝图像
%     figure
%     imshow(stick)
    
    skL = sketelon(stick); % 提取骨架
%     figure
%     imshow(skL)
    
    [ROWS,COLUMNS] = find(skL == 1);  % 提取骨架点坐标
    if isempty(ROWS)
        disp('something error...')
        return;
    else
        DDist = count_len(ROWS,COLUMNS);
        Pixel = sum(DDist);
        Lens(t) = K * Pixel + B;
    end
end