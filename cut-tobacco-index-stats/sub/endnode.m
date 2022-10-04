function [endpoint_index, nodepoint_index, rows, columns] = endnode(image)
%% Ѱ�ҽڵ�Ͷ˵����
%% input:
% image: Ŀ��ͼ��
%% output:
% endpoint_index:   �˵����
% nodepoint_index�� �ڵ����
% rows:    Ŀ���������ͼ��������,�����У�
% columns: Ŀ������� (ͼ������꣬������)
%% ����ʼ

%% ��ȡ����Ϊ1������
[rows, columns] = find(image == 1);
%% Ѱ�Ҷ˵�ͽڵ�
% �˵�
endpoint_index = [];    % �˵�����
Template_two = [0,1;1,2;2,3;3,4;4,5;5,6;6,7;0,7];  % ������������ģ��
Template_three = [0,1,7;1,2,3;3,4,5;5,6,7];        % ������������ģ��
% �ڵ�
nodepoint_index = [];   % �ڵ�����
for i = 1:length(rows)
    goal = [rows(i),columns(i)]; % Ŀ���
    [DomainValue, chain, field] = chaincode(image, goal); % ��������
    if length(chain) == 1 
        endpoint_index = [endpoint_index,i];   % �˵� 
    elseif length(chain) == 2 && ismember(chain,Template_two,'rows')
        endpoint_index = [endpoint_index,i];    % �˵�
    elseif length(chain) == 3 && ismember(chain,Template_three,'rows')
        endpoint_index = [endpoint_index,i];    % �˵�
    elseif length(chain) > 2
        nodepoint_index = [nodepoint_index,i];   % �ڵ�����
    end
end