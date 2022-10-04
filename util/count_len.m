function dist = count_len(rows,columns)
%% input
% rows: 纵坐标
% columns: 横坐标
%% output
% dist：像素值
%% 
x = [rows,columns];  % 组合横纵坐标
m_min = find(x(:,1)==min(x(:,1)));  % 纵坐标最小值
m_min = m_min(1); % 纵坐标最小值第一个点
xA = x(m_min,:); 
DD = x;
DD(m_min,:) = [];
xB = DD;
XU = 1:1:size(x,1);
XU(m_min) = [];
XUHAO = XU;
CC = zeros(size(x,1)-1);
dist = zeros(1, size(x,1)-1);
for k = 1:size(x,1)-1
    M = pdist2(xA,xB);
    [~,nn] = find(M==min(M));
    CC(k) = XUHAO(nn(1));
    XUHAO(nn(1)) = [];
    xA = xB(nn(1),:);
    xB(nn(1),:) = [];
    dist(k) = min(M);
   if min(M) > 3*mean(dist(dist > 0))
       dist(k) = 0;
   end
end