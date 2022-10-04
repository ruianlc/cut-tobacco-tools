clc;clear;close all

rootdir = 'G:\workspace\chaungheyi\02. 业务\01. 烟草项目\01. 原料_制丝_卷包\00. 数据源\卷包制丝\叶片结构-烟丝结构12月份分析\综合测试台标准样\';

Files = dir(fullfile(rootdir,'*.bmp'));
LengthFiles = length(Files);

TotalLen = cell(LengthFiles,1);

for kk = 1:LengthFiles
    
    source = imread(fullfile(rootdir,Files(kk).name)); % 读取每一张图片
    %source = I(200:end,:,:);
    
    rgb = reshape(source,size(source,1)*size(source,2),3);
    rgb(rgb(:,3)>200,:) = 0;
    rgb(rgb(:,3)~= 0,:) = 1;
    img = reshape(rgb,size(source,1),size(source,2),3);
    bw0 = img(:,:,1);
    bw = logical(bw0);
    bw2 = imfill(bw,'holes');
    adjustImage = bwareaopen(bw2, 500); % 剔除小面积图像
    
    regions = regionprops(adjustImage); % 计算连通区域，region = [{area,centroid,boundingbox}]
    [label,region_num] = bwlabel(adjustImage);
    
%     figure
%     imshow(adjustImage)
%     hold on
    
    singleImgLWA = zeros(region_num,3);
    tarea = cell(1,region_num); % 每个片烟面积
    
    for ss = 1:region_num
        bw = label;
        pos = regions(ss).BoundingBox;
        
        r1 = round(pos(2));
        c1 = round(pos(1));
        w = pos(3);
        h = pos(4);
        r2 = r1+h-1;
        c2 = c1+w-1;
        bw(bw ~= ss) = 0;  % 其他像素置0
        bw(bw == ss) = 1;  % 将等于i的像素置1
        subImg = bw(r1:r2,c1:c2);  % 每个片烟的像素区域
        
        A = regions(ss).Area;
        
%         figure
%         imshow(subImg)
        
        [r,c]=find(subImg==1);
        % 'a'是按面积算的最小矩形，如果按边长用'p'
        [rectx,recty,area,perimeter] = minboundrect(c,r,'p');
        
        len_1 = sqrt((rectx(1) - rectx(2))^2 + (recty(1) - recty(2))^2);
        len_2 = sqrt((rectx(2) - rectx(3))^2 + (recty(2) - recty(3))^2);
        L = max([len_1,len_2]); % 长度
        W = min([len_1,len_2]); % 宽度
        
        singleImgLWA(ss,:) = [L,W,A];
%         plot(rectx(1),recty(1),'g*')
%         plot(rectx(2),recty(2),'r*')
%         plot(rectx(3),recty(3),'m*')
%         
%         line(rectx,recty,'Color','r');
    end
    TotalLen{kk} = singleImgLWA;
        
end

LenResult = cell2mat(TotalLen);
writematrix(LenResult,'LWAresult.xlsx');