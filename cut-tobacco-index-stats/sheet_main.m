clc;clear;close all
%%
% --- ˵�� ---
% ��Ƭ��ͼ�����ģ����˿����ȡͼ����ÿ��Ƭ�̷ָ��������˿����
% 
% Programmer: Robin An, 2021.07.27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% TODO 
% ��С��Ӿء��������Բ����С����Բ����ͨ�����и�����

pathfile = 'G:\workspace\���ݴ洢\������˿���\������˿ʵ��\A\1.��\2021-09-08_16_29_15_082.bmp';
%pathfile = 'G:\workspace\���ݴ洢\������˿���\������˿ʵ��\A\��\2021-09-08_16_39_38_873.bmp';
%pathfile = 'G:\workspace\���ݴ洢\������˿���\������˿ʵ��\A\��\2021-09-08_17_30_17_381.bmp';
%pathfile = 'G:\workspace\���ݴ洢\������˿���\������˿ʵ��\A\С\2021-09-08_18_54_39_543.bmp';
%pathfile = 'G:\workspace\���ݴ洢\������˿���\������˿ʵ��\A\��\2021-09-09_10_08_35_786.bmp';

I = imread(pathfile);
%I = source(120:end,:,:);
%imtool(I)

% ��ֵͼ��
% rgb = reshape(I,size(I,1)*size(I,2),3);
% rgb(rgb(:,3)>90,:) = 0; % ����Bֵ����ȥ������
% rgb(rgb(:,3)~= 0,:) = 1;
% img = reshape(rgb,size(I,1),size(I,2),3);
% bw = img(:,:,1);
% BW = logical(bw);
% adjustImage = bwareaopen(BW, 600); % �޳�С���ͼ��

% nhood = true(9);  % �˲�����
% S = stdfilt(I,nhood);  % ��׼���˲�
% gray = mat2gray(S);  % �ҶȻ�
umb = rgb2gray(I); % ͼ��ĻҶȴ���
level = graythresh(umb);
bw = imbinarize(umb,level);  % ͼ��Ķ�ֵ������
bw1 = ~bw; % ���ÿ��Ƭ��
bw2 = imfill(bw1,'holes');
adjustImage = bwareaopen(bw2, 600); % �޳�С���ͼ��
% figure
% imshow(adjustImage)

%% ��ȡ����߽�
[boundaries,L] = bwboundaries(adjustImage);
% tcenter = cell(1,length(boundaries));
% tradius = zeros(1,length(boundaries));
area = zeros(1,length(boundaries));  % ����Ƭ�����
tarea = zeros(1,length(boundaries)); % ����Ƭ�����Բ���
rrate = zeros(1,length(boundaries)); % ����Ƭ��Բ����

% imshow(label2rgb(L, @jet, [.5 .5 .5]))
% hold on
% for k = 1:length(boundaries)
%     boundary = boundaries{k};
%     plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
% end

%% ��ȡÿ����ͨ����
regions = regionprops(adjustImage); % ������ͨ����region = [{area,centroid,boundingbox}]
[label,sheet_num] = bwlabel(adjustImage);

% �����ע��ͨ����
% tpos=regionprops(label,'BoundingBox');
% tcentroid = regionprops(label,'Centroid');
% for i = 1:sheet_num % ÿ��Ƭ��
%     rectangle('position',tpos(i).BoundingBox,'edgecolor','r', 'LineWidth',1.5);
%     line([tcentroid(i,1).Centroid(1,1), tcentroid(i,1).Centroid(1,1)],[tcentroid(i,1).Centroid(1,2)-35, tcentroid(i,1).Centroid(1,2)+35],'Color','b', 'LineWidth',2);
%     text(tcentroid(i,1).Centroid(1,1),tcentroid(i,1).Centroid(1,2), num2str(i),'Color', 'r');
% end

K = 0.1314;
B = 0;

PL = 8; % ������ʵ�ʳ��ȣ���˿���1mm����Ӧ����Ϊ12.1168����˿���2mm����Ӧ����Ϊ19

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
    subImg = bw(r1:r2,c1:c2);  % ÿ��Ƭ�̵���������
    area(i) = bwarea(subImg);  % ����Ƭ�����
    
    boundary = boundaries{i};
    [center,radius] = minboundcircle(boundary(:,2),boundary(:,1)); % ��С���Բ�뾶��Բ��
        
    % �������Բ
    % viscircles(center,radius);    
    % tcenter{i} = center;
    % tradius(i) = radius;
    tarea(i) = pi * (radius.^2);   
    rrate(i) = area(i) / tarea(i);
    
    %% Ƭ���и�
    % �Ե�����ͨ�������ģ����˿
    % �Կ��wΪ�ȷֶ�����10��������Ϊ3mm����Ϊ1����˿�Ŀ�ȣ�����ÿ����ͨ�����и������˿����
    stick_num = floor(w/PL);
    Stick = cell(1,stick_num); % ÿ����˿
    skL = cell(1,stick_num); % ÿ����˿�Ǽ�
    L = cell(1,stick_num); % һ��Ƭ�̵�������˿����
    
    %% Ƭ���и���˿����
    %savepic = strcat(cur_path,'\Figure\',num2str(i),'.png');
    figure
    %h0 = figure('color','b');
    %set(h0,'Visible','off');
    
    for j = 1:stick_num
        
        % ��ȡÿ����˿��Images{i}-->Stick{j}
        Stick{j} = subImg(:,(j-1)*PL+1:j*PL);
        
        % �жϴ˴�����ͨ����������и�������˿
        stick_regions = regionprops(Stick{j});
        stick_area = bwlabel(Stick{j});
        num_stick_regions = length(stick_regions);
        sStick = zeros(1,num_stick_regions); % �и���˿�ĳ��ȣ����ڶ�˿��
        
        if num_stick_regions == 1 % δ������˿
            skL{j} = sketelon(Stick{j}); % ��ȡ�Ǽ�
            [ROWS,COLUMNS] = find(skL{j} == 1);  % ��ȡ�Ǽܵ�����
            
            if size(ROWS,1) <= 20 % �޳�������С���
                continue;
            end
            DDist = count_len(ROWS,COLUMNS);
            DDist(DDist > 1.5) = [];
            Pixel = sum(DDist);
            L{j} = K * Pixel + B;
        else % ������˿
            
            for p = 1:num_stick_regions % ÿ����˿�����Ķ�˿����
                bw_stick = stick_area;
                pos_stick = stick_regions(j).BoundingBox; % ��Ӿ��ε�x��y�߽磨�½�����䣩��pos = [x,y,width,height]
                
                r1 = round(pos_stick(2));
                c1 = round(pos_stick(1));
                w = pos_stick(3);
                h = pos_stick(4);
                r2 = r1+h-1;
                c2 = c1+w-1;
                bw_stick(bw_stick ~= p) = 0;
                bw_stick(bw_stick == p) = 1;
                subStick{p} = bw_stick(r1:r2,c1:c2);
                
                %sskL{p} = bwmorph(subStick{p},'thin',inf);
                sketelon(subStick{p}); % ��ȡ�Ǽ�
                [ROWS,COLUMNS] = find(sskL{p} == 1);  % ��ȡ�Ǽܵ�����
                
                if size(ROWS,1) <= 20 % �޳�������С���
                    continue;
                end
                DDist = count_len(ROWS,COLUMNS);
                DDist(DDist > 1.5) = [];
                Pixel = sum(DDist);
                sStick(p) = K * Pixel + B;
            end
            L{j} = sStick(:);
        end
        
        subplot(1,stick_num,j);
        imshow(Stick{j});
    end
    %print(h0,'-dpng',savepic);
    LL{i} = [tarea(i),rrate(i),L{:}];
end

% [curpath,~] = fileparts(mfilename('fullpath')); % ��ǰĿ¼
% dsave = strcat(curpath,'\Ƭ���и���˿����ͳ��a.xlsx');
% writecell(LL,dsave,'Sheet',1,'Range','A1');
