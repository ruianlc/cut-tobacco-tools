clc;clear;close all
%%
% --- ˵�� ---
% ��Ƭ��ͼ�����ģ����˿����ȡͼ����ÿ��Ƭ�̷ָ��������˿����
% 
% Programmer: Robin An, 2021.07.27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TODO 
% ��С��Ӿء��������Բ����С����Բ����ͨ�����и�����

 pathfile = 'G:\workspace\���ݴ洢\��˿������ʵ��\Ƭ�̼��-2020.07.21\��Ƭ\20210721151138.bmp';
% pathfile = 'G:\workspace\���ݴ洢\��˿������ʵ��\Ƭ�̼��-2020.07.21\��Ƭ\20210721162841.bmp';
% pathfile = 'G:\workspace\���ݴ洢\��˿������ʵ��\20210721162841.bmp';
% pathfile = 'G:\workspace\���ݴ洢\��˿������ʵ��\Ƭ�̼��-2020.07.21\СƬ\20210721195251.bmp';
% pathfile = 'G:\workspace\���ݴ洢\��˿������ʵ��\Ƭ�̼��-2020.07.21\��Ƭ\20210721230318.bmp';

I = imread(pathfile);
figure
imshow(I)
title('ԭʼͼ��')
%I = source(120:end,:,:);
%imtool(I)

%% ��ֵͼ��
rgb = reshape(I,size(I,1)*size(I,2),3);
rgb(rgb(:,3)>150,:) = 0; % ȥ������
%rgb(rgb(:,1)<30,:) = 0;  % ȥ��Ӱ��Ԫ��

rgb(rgb(:,3)~= 0,:) = 1;
img = reshape(rgb,size(I,1),size(I,2),3);
bw = img(:,:,1);
BW = logical(bw);
figure
imshow(BW)
title('��ֵ��ͼ��')

% se = strel('disk',5); % �����Ǵ���һ���뾶Ϊ5��ƽ̹��Բ�̽ṹԪ��
% BW2 = imerode(BW1,se);% ��ʴ
% figure
% imshow(BW2)
% title('��ֵ����ʴͼ��')
% 
% montage({I,BW,BW1,BW2},'Size',[2 2])
adjustImage = bwareaopen(BW, 300); % �޳�С���ͼ��
figure
imshow(adjustImage)
title('ͼ���ֵ���޳�С���')

%% ��ȡ����߽�
[boundaries,L] = bwboundaries(adjustImage);
% tcenter = cell(1,length(boundaries));
% tradius = zeros(1,length(boundaries));
area = zeros(1,length(boundaries));  % ����Ƭ�����
tarea = zeros(1,length(boundaries)); % ����Ƭ�����Բ���
rrate = zeros(1,length(boundaries)); % ����Ƭ��Բ����

imshow(label2rgb(L, @jet, [.5 .5 .5]))
hold on
for k = 1:length(boundaries)
    boundary = boundaries{k};
    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
end

%% ��ȡÿ����ͨ����
regions = regionprops(adjustImage); % ������ͨ����region = [{area,centroid,boundingbox}]
[label,sheet_num] = bwlabel(adjustImage);

% �����ע��ͨ����
tpos=regionprops(label,'BoundingBox');
tcentroid = regionprops(label,'Centroid');
for i = 1:sheet_num % ÿ��Ƭ��
    rectangle('position',tpos(i).BoundingBox,'edgecolor','r', 'LineWidth',1.5);
    line([tcentroid(i,1).Centroid(1,1), tcentroid(i,1).Centroid(1,1)],[tcentroid(i,1).Centroid(1,2)-35, tcentroid(i,1).Centroid(1,2)+35],'Color','b', 'LineWidth',2);
    text(tcentroid(i,1).Centroid(1,1),tcentroid(i,1).Centroid(1,2), num2str(i),'Color', 'r');
end

K = 0.1314;
B = 0;
% K = 0.1336;
% B = -0.6188;

PL = 13; % ������ʵ�ʳ��ȣ���˿���1mm����Ӧ����Ϊ12.1168����˿���2mm����Ӧ����Ϊ19
Images = cell(1,sheet_num); % ÿ��ͼ������Ƭ����������

LL = cell(sheet_num,1);% ÿ��ͼ��Ƭ��ָ�꣺�����Բ���ʡ���˿����

[cur_path,~] = fileparts(mfilename('fullpath'));

for i = 1:sheet_num % ��ÿ��Ƭ�̷ֱ����
    bw = label;
    pos = regions(i).BoundingBox;
    
    r1 = round(pos(2));
    c1 = round(pos(1));
    w = pos(3);
    h = pos(4);
    r2 = r1+h-1;
    c2 = c1+w-1;
    bw(bw ~= i) = 0;  % ����������0
    bw(bw == i) = 1;  % ������i��������1
    Images{i} = bw(r1:r2,c1:c2);  % ÿ��Ƭ�̵���������
    nImage = ~Images{i};
        
    %% Ƭ���и�
    stick_num = floor(w/PL);
    
    %% Ƭ���и���˿����
    savepic = strcat(cur_path,'\Figure\',num2str(i),'.png');
    h0 = figure('color','w');
    set(h0,'Visible','off');
    
    for j = 1:stick_num-1 
        nImage(:,j*PL+1:j*PL+2) = 1;
    end
    imshow(nImage);
    print(h0,'-dpng',savepic);
end
