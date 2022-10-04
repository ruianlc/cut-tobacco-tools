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

rootdir = 'D:\workspace\数据存储\上烟卷接制丝\B模块切丝实验\在线制丝实验\A';
%rootdir = 'D:\workspace\数据存储\上烟卷接制丝\B模块切丝实验\在线制丝实验\B';
%rootdir = 'D:\workspace\数据存储\上烟卷接制丝\B模块切丝实验\在线制丝实验\B小样';

subdir = dir(rootdir);

% 大中小碎分别设定不同阈值
big_threshold = 500;
middle_threshold = 200;
small_threshold = 100;
detritus_threshold = 30;

% 面积
% K = 0.1051;
% B = 37.761;

% 长度
K = 0.1314;
B = 0;

% 6像素，0.8mm；8像素，1.0mm；10像素，1.2mm
PL = 8;%7.6104; % 像素与实际长度，切丝宽度1mm，对应像素为7.6104，以8个像素视为一根烟丝的宽度
% Images = cell(1,sheet_num); % 每个图像所有片烟所在区域

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
        source = imread(imgfile_path);
        I = source;
        rgb = reshape(I,size(I,1)*size(I,2),3);
        rgb(rgb(:,3)>90,:) = 0; % 根据B值差异去除背景
        rgb(rgb(:,3)~= 0,:) = 1;
        img = reshape(rgb,size(I,1),size(I,2),3);
        bw = img(:,:,1);
        BW = logical(bw);  
        adjustImage = bwareaopen(BW, threshold); % 剔除小面积图像

%         nhood = true(9);  % 滤波算子
%         S = stdfilt(I,nhood);  % 标准差滤波
%         gray = mat2gray(S);  % 灰度化
%         umb = rgb2gray(gray); % 图像的灰度处理
%         level = graythresh(umb);
%         bw = imbinarize(umb,level);  % 图像的二值化处理
%         bw1 = ~bw; % 获得每个片烟
%         bw2 = imfill(bw1,'holes');
%         adjustImage = bwareaopen(bw2, threshold); % 剔除小面积图像
%         figure
%         imshow(adjustImage)
        
        %% 获取区域边界
        [boundaries,L] = bwboundaries(adjustImage);
        area = zeros(1,length(boundaries));  % 所有片烟面积
        tarea = zeros(1,length(boundaries)); % 所有片烟外接圆面积
        rrate = zeros(1,length(boundaries)); % 所有片烟圆度率
   
        %% 提取每个连通区域
        regions = regionprops(adjustImage); % 计算连通区域，region = [{area,centroid,boundingbox}]
        [label,sheet_num] = bwlabel(adjustImage);

        Limage = cell(sheet_num,1);% 每张图像片烟指标：面积、圆度率、烟丝长度
        
        for s = 1:sheet_num % 对每个片烟分别操作
            bw = label;
            pos = regions(s).BoundingBox;
            
            r1 = round(pos(2));
            c1 = round(pos(1));
            w = pos(3);
            h = pos(4);
            r2 = r1+h-1;
            c2 = c1+w-1;
            bw(bw ~= s) = 0;  % 其他像素置0
            bw(bw == s) = 1;  % 将等于i的像素置1
            subImg = bw(r1:r2,c1:c2);  % 每个片烟的像素区域
            area(s) = bwarea(subImg);  % 计算片烟面积
            
%             figure
%             imshow(Images{s})
            
            boundary = boundaries{s};
            [center,radius] = minboundcircle(boundary(:,2),boundary(:,1)); % 最小外接圆半径和圆心
            
            % 画出外接圆
            % viscircles(center,radius);
            % tcenter{i} = center;
            % tradius(i) = radius;
            % tarea(s) = pi * (radius.^2);
            % rrate(s) = area(s) / tarea(s);
            
            %% 片烟切割
            % 对单个连通区域进行模拟切丝
            % 以宽度w为等分对象，以10个像素设为3mm宽，设为1根烟丝的宽度，计算每个连通区域切割出的烟丝数量
            stick_num = floor(w/PL);
            Lsheet = cell(1,stick_num); % 一个片烟的所有烟丝长度

            %% 片烟切割烟丝长度
            for t = 1:stick_num

                subStick = subImg(:,(t-1)*PL+1:t*PL);
                
                % 判断此处的连通区域个数，切割后产生碎丝
                stick_regions = regionprops(subStick);
                stick_area = bwlabel(subStick);
                num_stick_regions = length(stick_regions);
                sStick = zeros(1,num_stick_regions); % 切割烟丝的长度（存在断丝）
                
                if num_stick_regions == 1 % 未产生断丝
                    skL = sketelon(subStick); % 提取骨架
                    [ROWS,COLUMNS] = find(skL == 1);  % 提取骨架点坐标
                    
                    if size(ROWS,1) <= 20 % 剔除产生的小碎点
                        continue;
                    end
                    DDist = count_len(ROWS,COLUMNS);
                    DDist(DDist > 1.5) = []; % 剔除与边缘重合烟丝的异常值
                    Pixel = sum(DDist);
                    Lsheet{t} = K * Pixel + B;
                else % 产生断丝
                    
                    for p = 1:num_stick_regions % 每根切丝包含的断丝数量
                        bw_stick = stick_area;
                        pos_stick = stick_regions(t).BoundingBox; % 外接矩形的x，y边界（下界和区间），pos = [x,y,width,height]
                        
                        r1 = round(pos_stick(2));
                        c1 = round(pos_stick(1));
                        w = pos_stick(3);
                        h = pos_stick(4);
                        r2 = r1+h-1;
                        c2 = c1+w-1;
                        bw_stick(bw_stick ~= p) = 0;
                        bw_stick(bw_stick == p) = 1;
                        subsubStick = bw_stick(r1:r2,c1:c2);
                        
                        sskL = sketelon(subsubStick); % 提取骨架
                        [ROWS,COLUMNS] = find(sskL == 1);  % 提取骨架点坐标
                        
                        if size(ROWS,1) <= 20 % 剔除产生的小碎点
                            continue;
                        end
                        DDist = count_len(ROWS,COLUMNS);
                        DDist(DDist > 1.5) = []; % 剔除与边缘重合烟丝的异常值
                        Pixel = sum(DDist);
                        sStick(p) = K * Pixel + B;
                    end
                    Lsheet{t} = sStick(:);
                end
            end
            Limage{s} = [area(s),Lsheet{:}];
        end
        Ltotal = [Ltotal;Limage];
    end % j - 每个目录下的图像
    
    %% 统计结果：片烟面积、外接圆面积、圆度率
    [curpath,~] = fileparts(mfilename('fullpath')); % 当前目录
    dsave = strcat(curpath,'\results\片烟模拟切丝结果\','A.xlsx');
    
    stitle = {'序号','片烟面积（像素）','烟丝长度（mm）'};
    xuhao = 1:length(Ltotal);
    writecell(stitle,dsave,'Sheet',i-2,'Range','A1'); % 写入表头
    writematrix(xuhao',dsave,'Sheet',i-2,'Range','A2'); % 写入序号
    writecell(Ltotal,dsave,'Sheet',i-2,'Range','B2'); % 写入片烟面积
end %i - 子目录