clc;clear;close all
%%
% --- ˵�� ---
% ��ȡ�����Ŀ¼�µ�Ƭ��ͼ��
% ��ÿ��ͼ���е�Ƭ�̽���ģ����˿����ȡͼ����ÿ��Ƭ�̷ָ��������˿����
% 
% Programmer: Robin An, 2021.07.30
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO( DONE now )
% 1.��Ŀ¼��ȡ
% 2.count_len�����Ż�
% 3.�и������˿
% 4.����䣬Ƭ���ڲ����ڿն�����һ�����Ǵ˵㣬��Ҫ�ж�ÿ��Ƭ���ڲ��Ƿ���ڱպ����򣬼�Ƭ�����ڵıպ�������䣬��Ƭ���ڲ����ڵıպ����������
% 5.�޳��쳣ԭʼ����
% 6.��С���Բ
% 7.��С��Ӿ�

%rootdir = 'D:\workspace\datasource\������˿���\������˿ʵ��\A';
% rootdir = 'D:\workspace\datasource\������˿���\������˿ʵ��\B';
 rootdir = 'D:\workspace\datasource\������˿���\������˿ʵ��\BС��';

subdir = dir(rootdir);

% ����С��ֱ��趨��ͬ��ֵ
big_threshold = 800;
middle_threshold = 500;
small_threshold = 200;
detritus_threshold = 100;

for i = 1:length(subdir)
    
    % ����'.'��'..'Ŀ¼
    if(strcmp(subdir(i).name,'.')||...
            strcmp(subdir(i).name,'..')||...
            ~subdir(i).isdir())
        continue;
    end
    
    % ���ļ������Һ�׺Ϊbmp���ļ�
    subdirpath = fullfile( rootdir, subdir(i).name);
    subdirfiles = dir(strcat(subdirpath,'\','*.bmp'));
    num_img = length(subdirfiles); % ÿ�����ļ���ͼ������
    Ltotal = {};%cell(num_img,1); % cell�洢ÿ��ͼƬ��˿�����
    
    if(isequal(num_img,0))
        return;
    end
    
    if contains(subdir(i).name,'��')
        threshold = big_threshold;
    elseif contains(subdir(i).name,'��')
        threshold = middle_threshold;
    elseif contains(subdir(i).name,'С')
        threshold = small_threshold;
    elseif contains(subdir(i).name,'��')
        threshold = detritus_threshold;
    else
        threshold = middle_threshold;
    end
        
    for j = 1:num_img
        
        imgfile_path = fullfile(subdirpath, subdirfiles(j).name);
        I = imread(imgfile_path);
        
%         rgb = reshape(I,size(I,1)*size(I,2),3);
%         rgb(rgb(:,3)>110,:) = 0; % ����Bֵ����ȥ������
%         rgb(rgb(:,3)~= 0,:) = 1;
%         img = reshape(rgb,size(I,1),size(I,2),3);
%         bw = img(:,:,1);
%         BW = logical(bw);  
%         adjustImage = bwareaopen(BW, threshold); % �޳�С���ͼ��
        
        gray = rgb2gray(I); % ͼ��ĻҶȴ���
        level = graythresh(gray);
        bw = imbinarize(gray,level);  % ͼ��Ķ�ֵ������
        bw1 = ~bw; % ���ÿ��Ƭ��
        bw2 = imfill(bw1,'holes');
        adjustImage = bwareaopen(bw2, threshold); % �޳�С���ͼ��

        figure
        imshow(adjustImage)
    end % j - ÿ��Ŀ¼�µ�ͼ��
    
end %i - ��Ŀ¼