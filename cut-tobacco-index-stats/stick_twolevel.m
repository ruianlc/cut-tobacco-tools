clc;clear;close all
%%
% --- 说明 ---
% 读取多个文件夹下的叶丝图像，计算每张图像中的叶丝长度
% 
% Programmer: Robin An, 2021.07.26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rootdir = 'G:\workspace\数据存储\上烟制丝卷包\制丝段离线实验\叶丝检测-2021.07.22';
%rootdir = 'G:\workspace\数据存储\制丝段离线实验\少量图片测试';
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
        
        % 对单张图像的处理
        source = I(80:end-200,80:end-150,:);  % 图像切割
        % 图像滤波
        nhood = true(9);  % 滤波算子
        S = stdfilt(source,nhood);  % 标准差滤波
        gray = mat2gray(S);  % 灰度化
        umb = rgb2gray(gray); % 图像的灰度处理
        bw = imbinarize(umb,0.15);  % 图像的二值化处理，im2bw以rgb图像为输入，官方已不推荐使用；imbinarize以灰度图像为输入
        BW2 = imfill(bw,'holes');  % 填充图像
        adjustImage = bwareaopen(BW2, 100); % 剔除小面积图像
        contour = bwperim(adjustImage); % 提取轮廓
        regions = regionprops(adjustImage); % 计算联通区域
        area = bwlabel(adjustImage);
        % [area,Num] = bwlabel(adjustImage);
        
        
        % 提取每个烟丝图像
        % K=0.1007;
        % B=1.641;
        K = 1;
        B = 0;
        Num = length(regions);
        L = zeros(1,Num); % 单张图像中所有烟丝的长度
                
        for idx = 1:Num % idx为单张图像中每个烟丝的序号
            bw = area;
            pos = regions(idx).BoundingBox; % 最小外接矩形的x，y边界（下界和区间），pos = [x,y,width,height]
            r1 = round(pos(2));
            c1 = round(pos(1));
            w = pos(3);
            h = pos(4);
            r2 = r1+h-1;
            c2 = c1+w-1;
            bw(bw ~= idx) = 0;  % 其他像素置
            bw(bw == idx) = 1;  % 将等于i的像素置1
            Images{idx} = bw(r1:r2,c1:c2);  % 单个烟丝图像
            figure
            imshow(Images{idx})
            % Images{idx} = bw(r1:r2,c1:c2,:); 
            contour = bwperim(Images{idx});
            I1 = bwmorph(Images{idx},'skel',inf);
            [endpoint_index, nodepoint_index, rows, columns] = endnode(I1);
            [x1,y1] = find(contour == 1);
            [x2,y2] = find(I1 == 1);
            skL{idx} = sketelon(Images{idx}); % 提取骨架
            [ROWS,COLUMNS] = find(skL{idx} == 1);  % 提取骨架点坐标
            DDist = count_len(ROWS,COLUMNS);
            Pixel = sum(DDist);
            L(idx) = K * Pixel + B - 3;
        end
        
        % 保存结果
        Ntotal{j} = subdirfiles(j).name(1:end-4); % 每张图像名称，去掉文件后缀名
        Ltotal{j} = L;  % 每张图像包含的烟丝长度
    end
    [curpath,~] = fileparts(mfilename('fullpath')); % 当前目录
    dsave = strcat(curpath,'\烟丝长度统计.xlsx');
    writecell(Ntotal,dsave,'Sheet',i-2,'Range','A1'); % 第一列写入图像文件名
    writecell(Ltotal,dsave,'Sheet',i-2,'Range','B1'); % 写入烟丝长度
end