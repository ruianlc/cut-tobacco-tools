function [DomainValue, chain, field] = chaincode(image, goal)
%%  ����Ŀ��ǰ��Ŀ��ͼ���е�8����ֵ
%% inpput:
% image:   Ŀ��ͼ��
% goal:    Ŀ�����ͼ���е�����
%% output:
% DomainValue:   Ŀ�����ڵ�8����ֵ��1��8�У�
% chain:         Ŀ�����ڵ����뷽��
%% 8������
direction(1,:) = [0,1];    % 0 ����
direction(2,:) = [-1,1];   % 1 ����
direction(3,:) = [-1,0];   % 2 ����
direction(4,:) = [-1,-1];  % 3 ����
direction(5,:) = [0,-1];   % 4 ����
direction(6,:) = [1,-1];   % 5 ����
direction(7,:) = [1,0];    % 6 ����
direction(8,:) = [1,1];    % 7 ����
%% ����Ŀ���8����ֵ
DomainValue = zeros(1,8);  % ��ʼ��8����ֵ
field = zeros(8,2);
for i = 1:8
    field(i,:) = goal + direction(i,:);  % ����ÿ�������ͼ������
    coordinate = field(i,:);  % ����ÿ�������ͼ������
    % ���Ŀ������곬����Χ����
    if goal(1) > size(image,1)
        error('ͼ����������')
    elseif goal(2) > size(image,2)
        error('ͼ����������')
    end
    % �����������㳬��ͼ��Χ������Χ��ͼ������ֵ����0
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
