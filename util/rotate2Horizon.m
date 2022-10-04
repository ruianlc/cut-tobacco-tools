function [outBwImg] = rotate2Horizon(bwImage,segSize)
% --- 函数说明 ---
% 将任意角度的烟丝旋转至水平方向
%
% --- 输入 ---
% bwImage : 单根烟丝二值化图像
% segSize : 烟丝分段数目
%
% --- 输出 ---
% outBwImg : 旋转至水平方向的烟丝
%
% ------
% author: Robin An
% e-mail: ruianlc@gmail.com
% created: 2022-01-21,   using Matlab R2020b (9.9.0.1467703) 64-bit (win64)
% last modified by Robin An on 2022-1-21
% copyright 2022 - Robin An.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

skL = RemoveSmallBranches(bwImage);

padsize = 1;
padvalue = 0;

% 截取中间区域
if size(skL,1) > size(skL,2) % 行多
    startPos = round(size(skL,1)/segSize);
    endPos = round((segSize-1)*size(skL,1)/segSize);
    
    skL = skL(startPos:endPos,:);
    
else % 列多
    startPos = round(size(skL,2)/segSize);
    endPos = round((segSize-1)*size(skL,2)/segSize);
    
    skL = skL(:,startPos:endPos);
end

skLNew = padarray(skL, [padsize padsize], padvalue,'both'); % 消除当两端点均在边界上的异常
terminating_pts = find_skel_ends(skLNew);
%     figure,imshow(skL);
%     hold on; plot(terminating_pts(:,1),terminating_pts(:,2),'r*');

x_end_1 = terminating_pts(1,1); % 第一个端点的x坐标
y_end_1 = terminating_pts(1,2); % 第一个端点的y坐标
x_end_2 = terminating_pts(2,1); % 第二个端点的x坐标
y_end_2 = terminating_pts(2,2); % 第二个端点的y坐标
k = -(y_end_1 - y_end_2) / (x_end_1 - x_end_2);  % 斜率
angle = 360 * (atan(k)/(2*pi)); % 角度

% 将烟丝旋转至水平方向
if angle >=0
    % nSubImg = imrotate(subImg,-angle,'bilinear','crop'); % 与原始图像一样大小
    nSubImg = imrotate(bwImage,-angle);
else
    nSubImg = imrotate(bwImage,-180-angle);
end

outBwImg = nSubImg;

