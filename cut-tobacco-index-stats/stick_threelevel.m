clc;clear;close all

rootdir = 'G:\workspace\数据存储\上烟制丝卷包\在线实验-烟丝长度检测';

subdir = dir(rootdir);

% 长中短碎分别设定不同阈值
long_threshold = 300;
middle_threshold = 200;
short_threshold = 60;
detritus_threshold = 10;

for i = 1:length(subdir)
    
    % 跳过'.'和'..'目录
    if(strcmp(subdir(i).name,'.')||...
            strcmp(subdir(i).name,'..')||...
            ~subdir(i).isdir())
        continue;
    end
    
    subsubdir =dir(strcat(rootdir,'\',subdir(i).name));
    for j = 1:length(subsubdir)
        
        % 跳过'.'和'..'目录
        if(strcmp(subsubdir(j).name,'.')||...
                strcmp(subsubdir(j).name,'..')||...
                ~subsubdir(j).isdir())
            continue;
        end
        
        % 子文件夹下找后缀为bmp的文件
        subsubdirpath = fullfile( rootdir,subdir(i).name, subsubdir(j).name);
        subsubdirfiles = dir(strcat(subsubdirpath,'\','*.bmp'));
        num_img = length(subsubdirfiles); % 每个子文件夹图像总数
        Ltotal = cell(num_img,1); % cell存储每张图片烟丝的长度
        Ntotal = cell(num_img,1); % cell存储每张图片的名称
        if(isequal(num_img,0))
            return;
        end
        
        if contains(subsubdir(j).name,'长')
            threshold = long_threshold;
        elseif contains(subsubdir(j).name,'中')
            threshold = middle_threshold;
        elseif contains(subsubdir(j).name,'短')
            threshold = short_threshold;
        elseif contains(subsubdir(j).name,'碎')
            threshold = detritus_threshold;
        else
            threshold = middle_threshold;
        end
        
        for s = 1:num_img
            
            imgfile_path = fullfile(subsubdirpath, subsubdirfiles(s).name);
            I = imread(imgfile_path);
%             figure
%             imshow(I)
%             titleString = strcat('第',num2str(i),'个样本 -- ','文件夹：',subsubdir(j).name,' -- 第',num2str(s),'张图像');
%             title(titleString)
            
            source = I(1:end-100,200:end-200,:); %截去左下角阴影影响
%             figure
%             imshow(source)
%             title('处理图像')

            rgb = reshape(source,size(source,1)*size(source,2),3);
            rgb(rgb(:,3)>240,:) = 0; % 去除背景：B值大于240
            rgb(rgb(:,3)~= 0,:) = 1;
            img = reshape(rgb,size(source,1),size(source,2),3);
            bw0 = img(:,:,1);
            bw = logical(bw0);
%             figure
%             imshow(bw)
%             title('二值化图像')
            
            BW2 = imfill(bw,'holes');  % 填充图像
%             figure
%             imshow(BW2)
%             title('二值化填充图像')
            
            adjustImage = bwareaopen(BW2, threshold); % 剔除小面积图像
%             figure
%             imshow(adjustImage)
%             title('二值化剔除小面积图像')
%             hold on
            
            regions = regionprops(adjustImage); % 计算联通区域
            [labels,Num] = bwlabel(adjustImage);
            
            % 提取根个烟丝图像
%             K = 0.3262;
%             B = 0.7791;
            K = 0.1549; % 长度像素毫米转换
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
        
                skL = sketelon(stick); % 提取骨架
%                 figure
%                 imshow(skL)
                
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
            
            % 保存结果
            Ntotal{s} = strcat(subsubdir(j).name,subsubdirfiles(s).name(end-6:end)); % 每张图像名称，去掉文件后缀名
            Ltotal{s} = Lens;  % 每张图像包含的烟丝长度
        end
        [curpath,~] = fileparts(mfilename('fullpath'));
        desPath = strcat(curpath,'\reults\',subdir(i).name,'.xlsx');
        writecell(Ntotal,desPath,'Sheet',j-2,'Range','A1'); % 写入烟丝长度
        writecell(Ltotal,desPath,'Sheet',j-2,'Range','B1'); % 写入烟丝长度
    end
end