clc;clear;close all
%%
% [filename,pathname,~] = uigetfile({'*.bmp';'*.jpg'},'ѡ��ͼ��');  % ѡ���ļ���
% pathfile = strcat(pathname,'\',filename);  % ���·��+�ļ���
% I = imread(char(pathfile));  % ��ȡͼƬ
%pathfile = 'G:\workspace\���ݴ洢\������˿���\����ʵ��-��˿���ȼ��\#01\1.��\2021-09-07_11_01_00_041.bmp';
%pathfile = 'G:\workspace\���ݴ洢\������˿���\����ʵ��-��˿���ȼ��\#01\1.��\2021-09-07_11_22_25_344.bmp';
%pathfile = 'G:\workspace\���ݴ洢\������˿���\����ʵ��-��˿���ȼ��\#02\1.��\2021-09-07_13_06_04_216.bmp';
%pathfile = 'G:\workspace\���ݴ洢\������˿���\����ʵ��-��˿���ȼ��\#01\2.��\2021-09-07_11_45_26_790.bmp';
%pathfile = 'G:\workspace\���ݴ洢\������˿���\����ʵ��-��˿���ȼ��\#04\2.��\2021-09-07_14_26_28_168.bmp';
%pathfile = 'G:\workspace\���ݴ洢\������˿���\����ʵ��-��˿���ȼ��\#06\3.��\2021-09-07_15_34_46_247.bmp';
%pathfile = 'G:\workspace\���ݴ洢\������˿���\����ʵ��-��˿���ȼ��\#01\3.��\2021-09-07_11_43_36_383.bmp';
%pathfile = 'G:\workspace\���ݴ洢\������˿���\����ʵ��-��˿���ȼ��\#01\4.��\2021-09-07_12_00_02_420.bmp';
%pathfile = 'G:\workspace\���ݴ洢\������˿���\����ʵ��-��˿���ȼ��\#01\4.��\2021-09-07_13_03_35_649.bmp';

%pathfile = 'G:\workspace\���ݴ洢\������˿���\����ʵ��-��˿���ȼ��\#06\1.��\2021-09-07_15_39_52_775.bmp';
%pathfile = 'D:\workspace\�Ƽ���Ŀ\�Ϻ��̲�\���-��˿��Ŀ\������˿���\ʵ�鷽��\��Ӷ�ʵ��\��֧�𿪰ڷ���˿_1֧12-15���� - ����.jpg';
pathfile = 'D:\workspace\�Ƽ���Ŀ\�Ϻ��̲�\���-��˿��Ŀ\������˿���\ʵ�鷽��\��Ӷ�ʵ��\��֧�𿪰ڷ���˿_1֧12-15���� - ���� (2).jpg';
threshold = 1;

I = imread(pathfile);
figure
imshow(I)
%             titleString = strcat('��',num2str(i),'������ -- ','�ļ��У�',subsubdir(j).name,' -- ��',num2str(s),'��ͼ��');
%             title(titleString)

%             K = 0.3262;
%             B = 0.7791;
%             K = 0.1314; % �����ߴ�
%             B = 0;
%             pixel = 120;
%             L = K * pixel + B
%
%             100 13.1
%             120 15.8
%             150 19.7
%             180 23.7
%             200 26.3

source = I;
% source = I(1:end-100,200:end-200,:); %��ȥ���½���ӰӰ��
% figure
% imshow(source)
% title('����ͼ��')

rgb = reshape(source,size(source,1)*size(source,2),3);
rgb(rgb(:,3) > 200,:) = 0;
rgb(rgb(:,3)~= 0,:) = 1;
img = reshape(rgb,size(source,1),size(source,2),3);
bw0 = img(:,:,1);
bw = logical(bw0);


figure
imshow(bw)
title('��ֵ��ͼ��')

BW2 = imfill(bw,'holes');  % ���ͼ��
figure
imshow(BW2)
title('��ֵ�����ͼ��')
regions0 = regionprops(BW2);

adjustImage = bwareaopen(BW2, threshold); % �޳�С���ͼ��
figure
imshow(adjustImage)
title('��ֵ���޳�С���ͼ��')

regions = regionprops(adjustImage); % ������ͨ����
[labels,Num] = bwlabel(adjustImage);

% tpos=regionprops(labels,'BoundingBox');
% tcentroid = regionprops(labels,'Centroid');
% for i = 1:Num % ÿ��Ƭ��
%     rectangle('position',tpos(i).BoundingBox,'edgecolor','r', 'LineWidth',1.5);
%     line([tcentroid(i,1).Centroid(1,1), tcentroid(i,1).Centroid(1,1)],[tcentroid(i,1).Centroid(1,2)-35, tcentroid(i,1).Centroid(1,2)+35],'Color','b', 'LineWidth',2);
%     text(tcentroid(i,1).Centroid(1,1),tcentroid(i,1).Centroid(1,2), num2str(i),'Color', 'r');
% end

% ��ȡ������˿ͼ��
% K = 0.3262;
% B = 0.7791;
K = 1;
B = 0;
Lens = zeros(1,Num); % ����ͼ����������˿�ĳ���

for t = 1:Num % idxΪ����ͼ����ÿ����˿�����
    tlabel = labels;
    
    pos = regions(t).BoundingBox; % ��С��Ӿ��ε�x��y�߽磨�½�����䣩��pos = [x,y,width,height]
    r1 = round(pos(2));
    c1 = round(pos(1));
    w = pos(3);
    h = pos(4);
    r2 = r1+h-1;
    c2 = c1+w-1;
    tlabel(tlabel ~= t) = 0;  % ����������
    tlabel(tlabel == t) = 1;  % ������i��������1
    stick = tlabel(r1:r2,c1:c2);  % ������˿ͼ��
%     figure
%     imshow(stick)
    
    skL = sketelon(stick); % ��ȡ�Ǽ�
%     figure
%     imshow(skL)
    
    [ROWS,COLUMNS] = find(skL == 1);  % ��ȡ�Ǽܵ�����
    if isempty(ROWS)
        disp('something error...')
        return;
    else
        DDist = count_len(ROWS,COLUMNS);
        Pixel = sum(DDist);
        Lens(t) = K * Pixel + B;
    end
end