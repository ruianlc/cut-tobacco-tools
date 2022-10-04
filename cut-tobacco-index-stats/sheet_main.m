clc;clear;close all
%%
% --- 说明 ---
% 对片烟图像进行模拟切丝，获取图像中每个片烟分割出来的烟丝长度
% 
% Programmer: Robin An, 2021.07.27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TODO 
% 最小外接矩、最大内切圆、最小外切圆、连通区域切割、不填充

pathfile = 'G:\workspace\数据存储\上烟制丝卷包\在线制丝实验\A\1.大\2021-09-08_16_29_15_082.bmp';
%pathfile = 'G:\workspace\数据存储\上烟制丝卷包\在线制丝实验\A\大\2021-09-08_16_39_38_873.bmp';
%pathfile = 'G:\workspace\数据存储\上烟制丝卷包\在线制丝实验\A\中\2021-09-08_17_30_17_381.bmp';
%pathfile = 'G:\workspace\数据存储\上烟制丝卷包\在线制丝实验\A\小\2021-09-08_18_54_39_543.bmp';
%pathfile = 'G:\workspace\数据存储\上烟制丝卷包\在线制丝实验\A\碎\2021-09-09_10_08_35_786.bmp';

I = imread(pathfile);
%I = source(120:end,:,:);
%imtool(I)

% 二值图像
% rgb = reshape(I,size(I,1)*size(I,2),3);
% rgb(rgb(:,3)>90,:) = 0; % 根据B值差异去除背景
% rgb(rgb(:,3)~= 0,:) = 1;
% img = reshape(rgb,size(I,1),size(I,2),3);
% bw = img(:,:,1);
% BW = logical(bw);
% adjustImage = bwareaopen(BW, 600); % 剔除小面积图像

% nhood = true(9);  % 滤波算子
% S = stdfilt(I,nhood);  % 标准差滤波
% gray = mat2gray(S);  % 灰度化
umb = rgb2gray(I); % 图像的灰度处理
level = graythresh(umb);
bw = imbinarize(umb,level);  % 图像的二值化处理
bw1 = ~bw; % 获得每个片烟
bw2 = imfill(bw1,'holes');
adjustImage = bwareaopen(bw2, 600); % 剔除小面积图像
% figure
% imshow(adjustImage)

%% 获取区域边界
[boundaries,L] = bwboundaries(adjustImage);
% tcenter = cell(1,length(boundaries));
% tradius = zeros(1,length(boundaries));
area = zeros(1,length(boundaries));  % 所有片烟面积
tarea = zeros(1,length(boundaries)); % 所有片烟外接圆面积
rrate = zeros(1,length(boundaries)); % 所有片烟圆度率

% imshow(label2rgb(L, @jet, [.5 .5 .5]))
% hold on
% for k = 1:length(boundaries)
%     boundary = boundaries{k};
%     plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
% end

%% 提取每个连通区域
regions = regionprops(adjustImage); % 计算连通区域，region = [{area,centroid,boundingbox}]
[label,sheet_num] = bwlabel(adjustImage);

% 方框标注连通区域
% tpos=regionprops(label,'BoundingBox');
% tcentroid = regionprops(label,'Centroid');
% for i = 1:sheet_num % 每个片烟
%     rectangle('position',tpos(i).BoundingBox,'edgecolor','r', 'LineWidth',1.5);
%     line([tcentroid(i,1).Centroid(1,1), tcentroid(i,1).Centroid(1,1)],[tcentroid(i,1).Centroid(1,2)-35, tcentroid(i,1).Centroid(1,2)+35],'Color','b', 'LineWidth',2);
%     text(tcentroid(i,1).Centroid(1,1),tcentroid(i,1).Centroid(1,2), num2str(i),'Color', 'r');
% end

K = 0.1314;
B = 0;

PL = 8; % 像素与实际长度，切丝宽度1mm，对应像素为12.1168；切丝宽度2mm，对应像素为19

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
    subImg = bw(r1:r2,c1:c2);  % 每个片烟的像素区域
    area(i) = bwarea(subImg);  % 计算片烟面积
    
    boundary = boundaries{i};
    [center,radius] = minboundcircle(boundary(:,2),boundary(:,1)); % 最小外接圆半径和圆心
        
    % 画出外接圆
    % viscircles(center,radius);    
    % tcenter{i} = center;
    % tradius(i) = radius;
    tarea(i) = pi * (radius.^2);   
    rrate(i) = area(i) / tarea(i);
    
    %% 片烟切割
    % 对单个连通区域进行模拟切丝
    % 以宽度w为等分对象，以10个像素设为3mm宽，设为1根烟丝的宽度，计算每个连通区域切割出的烟丝数量
    stick_num = floor(w/PL);
    Stick = cell(1,stick_num); % 每根烟丝
    skL = cell(1,stick_num); % 每根烟丝骨架
    L = cell(1,stick_num); % 一个片烟的所有烟丝长度
    
    %% 片烟切割烟丝长度
    %savepic = strcat(cur_path,'\Figure\',num2str(i),'.png');
    figure
    %h0 = figure('color','b');
    %set(h0,'Visible','off');
    
    for j = 1:stick_num
        
        % 获取每根烟丝，Images{i}-->Stick{j}
        Stick{j} = subImg(:,(j-1)*PL+1:j*PL);
        
        % 判断此处的连通区域个数，切割后产生碎丝
        stick_regions = regionprops(Stick{j});
        stick_area = bwlabel(Stick{j});
        num_stick_regions = length(stick_regions);
        sStick = zeros(1,num_stick_regions); % 切割烟丝的长度（存在断丝）
        
        if num_stick_regions == 1 % 未产生断丝
            skL{j} = sketelon(Stick{j}); % 提取骨架
            [ROWS,COLUMNS] = find(skL{j} == 1);  % 提取骨架点坐标
            
            if size(ROWS,1) <= 20 % 剔除产生的小碎点
                continue;
            end
            DDist = count_len(ROWS,COLUMNS);
            DDist(DDist > 1.5) = [];
            Pixel = sum(DDist);
            L{j} = K * Pixel + B;
        else % 产生断丝
            
            for p = 1:num_stick_regions % 每根切丝包含的断丝数量
                bw_stick = stick_area;
                pos_stick = stick_regions(j).BoundingBox; % 外接矩形的x，y边界（下界和区间），pos = [x,y,width,height]
                
                r1 = round(pos_stick(2));
                c1 = round(pos_stick(1));
                w = pos_stick(3);
                h = pos_stick(4);
                r2 = r1+h-1;
                c2 = c1+w-1;
                bw_stick(bw_stick ~= p) = 0;
                bw_stick(bw_stick == p) = 1;
                subStick{p} = bw_stick(r1:r2,c1:c2);
                
                %sskL{p} = bwmorph(subStick{p},'thin',inf);
                sketelon(subStick{p}); % 提取骨架
                [ROWS,COLUMNS] = find(sskL{p} == 1);  % 提取骨架点坐标
                
                if size(ROWS,1) <= 20 % 剔除产生的小碎点
                    continue;
                end
                DDist = count_len(ROWS,COLUMNS);
                DDist(DDist > 1.5) = [];
                Pixel = sum(DDist);
                sStick(p) = K * Pixel + B;
            end
            L{j} = sStick(:);
        end
        
        subplot(1,stick_num,j);
        imshow(Stick{j});
    end
    %print(h0,'-dpng',savepic);
    LL{i} = [tarea(i),rrate(i),L{:}];
end

% [curpath,~] = fileparts(mfilename('fullpath')); % 当前目录
% dsave = strcat(curpath,'\片烟切割烟丝长度统计a.xlsx');
% writecell(LL,dsave,'Sheet',1,'Range','A1');
