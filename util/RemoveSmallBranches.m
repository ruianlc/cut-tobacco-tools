function [skelD] = RemoveSmallBranches(bw)
% --- ����˵�� ---
% ȥ���Ǽ�ͼ�еķֲ�֧
%
% --- ���� ---
% bw : ��ֵͼ
%
% --- ��� ---
% skelD : ������֧������ɹǼ�
%
% Programmer: Robin An, 2021-09-24
% last modified by Robin An on 2021-09-24
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%I got a better starting image with the 'thin' option than the 'skel' option

In = imfill(bw,'holes'); % Ѫ��ͼ���м䲿�����
I = bwmorph(In,'thin',Inf);

%Alternative splitting method to 'branchpoint'
%Use convolution to identify points with more than 2 neighboring pixels
% �зֲ棬����ζ�Ž�����г���2�����ϵ���������

% ����ˣ���Сһ����3*3��5*5��7*7
% ���ѧϰ��ʹ���ݶ��½��ͷ��򴫲��Զ�ѧϰ�������

filter = [1 1 1;
          1 0 1;
          1 1 1];

% returns the central part of the convolution that is the same size as A.
% C = conv2(A,B) ���ؾ��� A �� B �Ķ�ά���
% C = conv2(A,B,'same') ���ؾ���д�С�� A ��ͬ�����Ĳ���
convtmp = conv2(double(I), filter, 'same');
I_disconnect = I & ~(I & convtmp>2);

% CC = bwconncomp(BW)
% % ȷ��ͼ���е��������������������������������Ϊ 0����
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
% ɾ����ɢ����
skelD = bwmorph(I, 'spur');
% figure; imshow(skelD);

