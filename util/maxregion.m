function img = maxregion(I)
% --- 函数说明 ---
% 保留面积最大的区域
%
% 1.输入为二值化图，返回二值化图
% 2.输入为灰度图，返回灰度图
% 3.输入为rgb图，返回rgb图
%
% --- 输入 ---
% I : 图像
%
% --- 输出 ---
% img : 保留的最大区域图像
%
% ------
% author: Robin An
% e-mail: ruianlc@gmail.com
% created: 2021-09-23,   using Matlab R2020b (9.9.0.1467703) 64-bit (win64)
% last modified by Robin An on 2021-10-23
% copyright 2021 - Robin An.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lenI = length(size(I));

if lenI == 2
    if ~islogical(I)         % 灰度图
        imBw = imbinarize(I,greythresh(I));                       
    else                     % 二值化图
        imBw = I;
    end
    
    imLabel = bwlabel(imBw);
    stats = regionprops(imLabel,'Area');
    area = cat(1,stats.Area);
    iid = find(area == max(area));
    
    imLabel(imLabel ~= iid) = 0;
    imLabel(imLabel == iid) = 1;
    
    if ~islogical(I)
        img = imLabel.*I;
    else
        img = imLabel;
    end
else            % rgb图
    rgbImage = I;
    grayImage = min(rgbImage, [], 3); % 获取最大区域

    threshold = graythresh(grayImage);
    mask = imbinarize(grayImage,threshold);
    
    % 提取面积最大的区域
    mask = bwareafilt(~mask, 1);
    mask = imfill(mask, 'holes');

    redChannel = rgbImage(:, :, 1);
    greenChannel = rgbImage(:, :, 2);
    blueChannel = rgbImage(:, :, 3);
    
    % 三个通道的灰度化均值
    meanGLR = mean(rgbImage(redChannel > threshold));
    meanGLG = mean(rgbImage(greenChannel > threshold));
    meanGLB = mean(rgbImage(blueChannel > threshold));

    % mask以外的区域使用平均灰度填充
    redChannel(~mask) = meanGLR;
    greenChannel(~mask) = meanGLG;
    blueChannel(~mask) = meanGLB;
    
    img = cat(3, redChannel, greenChannel, blueChannel);

end
