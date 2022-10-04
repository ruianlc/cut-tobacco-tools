function [DomainValue, chain, field] = chaincode(image, goal)
%%  计算目标前在目标图像中的8领域值
%% inpput:
% image:   目标图像
% goal:    目标点在图像中的坐标
%% output:
% DomainValue:   目标点对于的8领域值（1行8列）
% chain:         目标点对于的链码方向
%% 8领域方向
direction(1,:) = [0,1];    % 0 方向
direction(2,:) = [-1,1];   % 1 方向
direction(3,:) = [-1,0];   % 2 方向
direction(4,:) = [-1,-1];  % 3 方向
direction(5,:) = [0,-1];   % 4 方向
direction(6,:) = [1,-1];   % 5 方向
direction(7,:) = [1,0];    % 6 方向
direction(8,:) = [1,1];    % 7 方向
%% 计算目标点8领域值
DomainValue = zeros(1,8);  % 初始化8领域值
field = zeros(8,2);
for i = 1:8
    field(i,:) = goal + direction(i,:);  % 计算每个领域点图像坐标
    coordinate = field(i,:);  % 计算每个领域点图像坐标
    % 如果目标点坐标超出范围报错
    if goal(1) > size(image,1)
        error('图像行数超出')
    elseif goal(2) > size(image,2)
        error('图像列数超出')
    end
    % 如果领域坐标点超出图像范围超出范围，图像领域值记作0
    if ((coordinate(1) <= 0) || (coordinate(2) <= 0) || ...
            (coordinate(1) > size(image,1)) || ...
            (coordinate(2) > size(image,2)))
        DomainValue(i) = 0;
    else
        DomainValue(i) = image(coordinate(1),coordinate(2));
    end
end
%% 
chain = find(DomainValue == 1) - 1;
