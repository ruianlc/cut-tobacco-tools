function [endpoint_index, nodepoint_index, rows, columns] = endnode(image)
%% 寻找节点和端点程序
%% input:
% image: 目标图像
%% output:
% endpoint_index:   端点序号
% nodepoint_index： 节点序号
% rows:    目标点行数（图像纵坐标,矩阵行）
% columns: 目标点列数 (图像横坐标，矩阵列)
%% 程序开始

%% 提取像素为1的坐标
[rows, columns] = find(image == 1);
%% 寻找端点和节点
% 端点
endpoint_index = [];    % 端点索引
Template_two = [0,1;1,2;2,3;3,4;4,5;5,6;6,7;0,7];  % 领域内两个点模板
Template_three = [0,1,7;1,2,3;3,4,5;5,6,7];        % 领域内三个点模板
% 节点
nodepoint_index = [];   % 节点索引
for i = 1:length(rows)
    goal = [rows(i),columns(i)]; % 目标点
    [DomainValue, chain, field] = chaincode(image, goal); % 计算链码
    if length(chain) == 1 
        endpoint_index = [endpoint_index,i];   % 端点 
    elseif length(chain) == 2 && ismember(chain,Template_two,'rows')
        endpoint_index = [endpoint_index,i];    % 端点
    elseif length(chain) == 3 && ismember(chain,Template_three,'rows')
        endpoint_index = [endpoint_index,i];    % 端点
    elseif length(chain) > 2
        nodepoint_index = [nodepoint_index,i];   % 节点索引
    end
end