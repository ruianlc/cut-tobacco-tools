function dist = count_len(rows,columns)
%% input
% rows: ������
% columns: ������
%% output
% dist������ֵ
%% 
x = [rows,columns];  % ��Ϻ�������
num_pixel = size(x,1);
num_dis = num_pixel - 1;
tmp_dist = zeros(1,num_dis);

m_min = find(x(:,1)==min(x(:,1)));  % ��������Сֵ
m_min = m_min(1); % ��������Сֵ��һ����
xA = x(m_min,:); % ��������Сֵ������
DD = x;
DD(m_min,:) = [];
xB = DD; % ����������Сֵ�����������������
% ZZ = m_min;
% XU = [1:1:size(x,1)];
XU = 1:num_pixel;
XU(m_min) = [];
XUHAO = XU;
for k = 1:num_dis
    M = pdist2(xA,xB);
    [~,nn] = find(M==min(M));
    % CC(k) = XUHAO(nn(1));
    XUHAO(nn(1)) = [];
    xA = xB(nn(1),:);
    xB(nn(1),:) = [];
    tmp_dist(k) = min(M);
%    if min(M) > 3*mean(dist)
%        dist(k) = 0;
%    end
end
dist = tmp_dist;
% redun = 1;
%  CC = [ZZ,CC];
%  plot(m1,n1,'ro')
%  hold on
%  plot(m1(CC),n1(CC))