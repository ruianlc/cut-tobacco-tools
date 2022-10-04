function [outBwImg] = rotate2Horizon(bwImage,segSize)
% --- ����˵�� ---
% ������Ƕȵ���˿��ת��ˮƽ����
%
% --- ���� ---
% bwImage : ������˿��ֵ��ͼ��
% segSize : ��˿�ֶ���Ŀ
%
% --- ��� ---
% outBwImg : ��ת��ˮƽ�������˿
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

% ��ȡ�м�����
if size(skL,1) > size(skL,2) % �ж�
    startPos = round(size(skL,1)/segSize);
    endPos = round((segSize-1)*size(skL,1)/segSize);
    
    skL = skL(startPos:endPos,:);
    
else % �ж�
    startPos = round(size(skL,2)/segSize);
    endPos = round((segSize-1)*size(skL,2)/segSize);
    
    skL = skL(:,startPos:endPos);
end

skLNew = padarray(skL, [padsize padsize], padvalue,'both'); % ���������˵���ڱ߽��ϵ��쳣
terminating_pts = find_skel_ends(skLNew);
%     figure,imshow(skL);
%     hold on; plot(terminating_pts(:,1),terminating_pts(:,2),'r*');

x_end_1 = terminating_pts(1,1); % ��һ���˵��x����
y_end_1 = terminating_pts(1,2); % ��һ���˵��y����
x_end_2 = terminating_pts(2,1); % �ڶ����˵��x����
y_end_2 = terminating_pts(2,2); % �ڶ����˵��y����
k = -(y_end_1 - y_end_2) / (x_end_1 - x_end_2);  % б��
angle = 360 * (atan(k)/(2*pi)); % �Ƕ�

% ����˿��ת��ˮƽ����
if angle >=0
    % nSubImg = imrotate(subImg,-angle,'bilinear','crop'); % ��ԭʼͼ��һ����С
    nSubImg = imrotate(bwImage,-angle);
else
    nSubImg = imrotate(bwImage,-180-angle);
end

outBwImg = nSubImg;

