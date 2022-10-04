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

% rootdir = 'G:\workspace\数据存储\上烟制丝卷包\制丝段离线实验\片烟检测-2020.07.21';
% rootdir = 'G:\workspace\数据存储\上烟制丝卷包\制丝段离线实验\少量图片测试';
rootdir = 'G:\workspace\数据存储\上烟制丝卷包\在线制丝实验\A';
subdir = dir(rootdir);

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
    Ltotal = cell(num_img,1); % cell存储每张图片烟丝的长度
    Ntotal = cell(num_img,1); % cell存储每张图片的名称
    
    if(isequal(num_img,0))
        return;
    end
    
    for j = 1:num_img
        
        imgfile_path = fullfile(subdirpath, subdirfiles(j).name);
        I = imread(imgfile_path);
        
       %% 二值图像    
        gray = rgb2gray(I); % 图像的灰度处理
        level = graythresh(gray);
        bw = imbinarize(gray,level);  % 图像的二值化处理
        bw1 = ~bw; % 获得每个片烟
        adjustImage = bwareaopen(bw1, 300); % 剔除小面积图像
       
        % 提取每个连通区域
        regions = regionprops(adjustImage); % 计算联通区域，region = [{area,centroid,boundingbox}]
        [area,sheet_num] = bwlabel(adjustImage);
        tpos=regionprops(area,'BoundingBox');
        tcentroid = regionprops(area,'Centroid');
        
        % 方框标注连通区域
        % for i = 1:sheet_num % 每个片烟
        %     rectangle('position',tpos(i).BoundingBox,'edgecolor','r', 'LineWidth',1.5);
        %     line([tcentroid(i,1).Centroid(1,1), tcentroid(i,1).Centroid(1,1)],[tcentroid(i,1).Centroid(1,2)-35, tcentroid(i,1).Centroid(1,2)+35],'Color','b', 'LineWidth',2);
        %     text(tcentroid(i,1).Centroid(1,1),tcentroid(i,1).Centroid(1,2), num2str(i),'Color', 'r');
        % end
        
        K = 0.1314;
        B = 0;
        
        PL = 7.6104; % 像素与实际长度，切丝宽度1mm，对应像素为7.6104
        Images = cell(1,sheet_num);
        LL = cell(1,sheet_num);% 每张图像所有烟丝长度
       
        % 获得每个连通区域的外接矩，得到长度
        for k = 1:sheet_num % 对每个片烟分布操作
            bw = area;
            pos = regions(k).BoundingBox; % 外接矩形的x，y边界（下界和区间），pos = [x,y,width,height]
            
            r1 = round(pos(2));
            c1 = round(pos(1));
            w = pos(3);
            h = pos(4);
            r2 = r1+h-1;
            c2 = c1+w-1;
            bw(bw ~= k) = 0;  % 其他像素置
            bw(bw == k) = 1;  % 将等于k的像素置1
            Images{k} = bw(r1:r2,c1:c2);  % 每个片烟的像素区域
            
            % 对单个连通区域进行模拟切丝。以宽度w为等分对象，以10个像素设为3mm宽，设为1根烟丝的宽度，计算每个连通区域切割出的烟丝数量
            stick_num = floor(w/PL);
            Stick = cell(1,stick_num); % 每根烟丝
            skL = cell(1,stick_num); % 每根烟丝骨架
            L = cell(1,stick_num); % 一个片烟的所有烟丝长度
            
            for t = 1:stick_num
                
                % 获取每根烟丝，Images{k}-->Stick{t}
                Stick{t} = Images{k}(:,(t-1)*PL+1:t*PL);
                
              %% TODO 判断此处的连通区域个数，切割后产生碎丝
                stick_regions = regionprops(Stick{t}); 
                stick_area = bwlabel(Stick{t});
                num_stick_regions = length(stick_regions);
                sStick = zeros(1,num_stick_regions);
                
                if num_stick_regions == 1 % 未产生断丝
                    skL{t} = sketelon(Stick{t}); % 提取骨架
                    [ROWS,COLUMNS] = find(skL{t} == 1);  % 提取骨架点坐标
                    
                    if size(ROWS,1) <= 20 % 剔除产生的小碎点
                        continue;
                    end
                    DDist = count_len(ROWS,COLUMNS);
                    Pixel = sum(DDist);
                    L{t} = K * Pixel + B;
                else % 产生断丝
                    
                    for p = 1:num_stick_regions % 每根切丝包含的断丝数量
                        bw_stick = stick_area;
                        pos_stick = stick_regions(k).BoundingBox; % 外接矩形的x，y边界（下界和区间），pos = [x,y,width,height]
                        
                        r1 = round(pos_stick(2));
                        c1 = round(pos_stick(1));
                        w = pos_stick(3);
                        h = pos_stick(4);
                        r2 = r1+h-1;
                        c2 = c1+w-1;
                        bw_stick(bw_stick ~= p) = 0;  
                        bw_stick(bw_stick == p) = 1; 
                        subStick{p} = bw_stick(r1:r2,c1:c2);  
                        
                        sskL{p} = sketelon(subStick{p}); % 提取骨架
                        [ROWS,COLUMNS] = find(sskL{p} == 1);  % 提取骨架点坐标
                        
                        if size(ROWS,1) <= 20 % 剔除产生的小碎点
                            continue;
                        end
                        DDist = count_len(ROWS,COLUMNS);
                        Pixel = sum(DDist);
                        sStick(p) = K * Pixel + B;
                    end
                    L{t} = sStick(1,:);
                end
            end % t - 每个片烟的烟丝
            LL{k} = [L{:}];
        end %k - 每张图像中的片烟 
        Ntotal{j} = subdirfiles(j).name(1:end-4); 
        Ltotal{j} = [LL{:}];
    end % j - 每个目录下的图像
    
    [curpath,~] = fileparts(mfilename('fullpath')); % 当前目录
    dsave = strcat(curpath,'\片烟切割烟丝长度统计.xlsx');
    writecell(Ntotal,dsave,'Sheet',i-2,'Range','A1'); % 第一列写入图像文件名
    writecell(Ltotal,dsave,'Sheet',i-2,'Range','B1'); % 写入烟丝长度
end %i - 子目录