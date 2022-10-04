function [skelD] = RemoveSmallBranches(bw)
% --- 函数说明 ---
% 去除骨架图中的分叉支
%
% --- 输入 ---
% bw : 二值图
%
% --- 输出 ---
% skelD : 消除分支后的主干骨架
%
% Programmer: Robin An, 2021-09-24
% last modified by Robin An on 2021-09-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%I got a better starting image with the 'thin' option than the 'skel' option

In = imfill(bw,'holes'); % 血管图像中间部分填充
I = bwmorph(In,'thin',Inf);

%Alternative splitting method to 'branchpoint'
%Use convolution to identify points with more than 2 neighboring pixels
% 有分叉，则意味着交叉点有超过2个以上的相邻像素

% 卷积核，大小一般是3*3、5*5、7*7
% 深度学习中使用梯度下降和反向传播自动学习出卷积核

filter = [1 1 1;
          1 0 1;
          1 1 1];

% returns the central part of the convolution that is the same size as A.
% C = conv2(A,B) 返回矩阵 A 和 B 的二维卷积
% C = conv2(A,B,'same') 返回卷积中大小与 A 相同的中心部分
convtmp = conv2(double(I), filter, 'same');
I_disconnect = I & ~(I & convtmp>2);

% CC = bwconncomp(BW)
% % 确定图像中的最大分量并将其擦除（将所有像素设置为 0）。
% numPixels = cellfun(@numel,CC.PixelIdxList);
% [biggest,idx] = max(numPixels);
% BW(CC.PixelIdxList{idx}) = 0;

cc = bwconncomp(I_disconnect);
numPixels = cellfun(@numel,cc.PixelIdxList);
[sorted_px, ind] = sort(numPixels);

%Remove components shorter than threshold
threshold  = 10;
for ii=ind(sorted_px<threshold)
    cur_comp = cc.PixelIdxList{ii};
    I(cur_comp) = 0; 

    %Before removing component, check whether image is still connected
    full_cc = bwconncomp(I);
    if full_cc.NumObjects>1
        I(cur_comp) = 1; 
    end
end

%Clean up left over spurs
% 删除杂散像素
skelD = bwmorph(I, 'spur');
% figure; imshow(skelD);

