function skL = sketelon(bws)
%% 骨架提取程序
%% input:
% bws：单个区域片烟二值图像（矩形切割后）
%% output:
% skL：最终提取的骨架
%% 程序开始
sk = bwmorph(bws,'skel',Inf);  % 提取骨架
bm = ones(3); % 生成3*3的全为1的矩阵
bm(2,2) = 0;  % 将中间的数值设置为0
[ro,la] = size(bws); % 二值图像的大小
ze1 = zeros(ro,1);  % 生成同二值图相同行数的全0列向量
ze2 = zeros(1,la+2); % 生成一行二值图列数+2的行向量
sk1 = [ze1,sk,ze1];  % 给骨架图像左右两边各加一列0向量
sk2 = [ze2;sk1;ze2];  % 给骨架上下各加一行0向量
sumsk = zeros(ro,la);  % 生成全0的同二值图像大小的0矩阵
%% 寻找节点
for i = 2:ro
    for j = 2:la
        SK = sk2(i-1:i+1,j-1:j+1);  % 每次提取3*3的矩阵
        SK1 = SK .* bm;  % 将该矩阵中心数值置0
        su = sum(sum(SK1));  % 求该矩阵的和
        if (su > 2) && (sk2(i,j) ~= 0)  % 如果和大于2并且该点值不等于0
            sumsk(i-1,j-1) = su;  % 节点二值化图
        end
    end
end
[front1,fn1] = bwlabel(sumsk); % 图像标记
skk = sk;  % 骨架图像
szo = zeros(ro,la);  % 全零二值图像
for m = 1:fn1
    sk0 = skk;  % 骨架图
    sk00 = skk;  % 骨架图
    sk0(find(front1==m)) = 0;  % 将骨架图像中该节点删除
    [front2,fn2] = bwlabel(sk0); % 标记删除该节点后的图像
    if fn2 < 3  % 如果分支小于3个，程序跳出
        continue
    end
    farea = front2(:);  % 将矩阵转成向量
    [fcounts,~] = histc(farea,1:fn2);  % 统计各分支的长度
    [~,tn] = sort(fcounts); % 对长度进行排序（由小到大）
    %[~,st]=min(fcounts);
    for i = 1:length(tn)-2
        sk00(find(front2 == tn(i))) = 0;  % 留下最长的两条分支
    end
    ss = skk - sk00;  % 分支图像
    szo = ss + szo;  % 每个节点分支
end
szo(find(szo ~= 0)) = 1;
skL = skk - szo;  % 去除分支后的节点图像
skL = logical(skL);
skL = bwmorph(skL,'spur',1);  % 再剔除小的分支
