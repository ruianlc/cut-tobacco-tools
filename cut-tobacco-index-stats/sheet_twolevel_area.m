clc;clear;close all
%%
% --- ˵�� ---
% ��ȡ�����Ŀ¼�µ�Ƭ��ͼ��
% ��ÿ��ͼ���е�Ƭ�̽���ģ����˿����ȡͼ����ÿ��Ƭ�̷ָ��������˿����
% 
% Programmer: Robin An, 2021.07.28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TODO
% DONE
% 1.��Ŀ¼��ȡ
% 2.count_len�����Ż�
% 3.�и������˿
% 4.����䣬Ƭ���ڲ����ڿն�����һ�����Ǵ˵㣬��Ҫ�ж�ÿ��Ƭ���ڲ��Ƿ���ڱպ����򣬼�Ƭ�����ڵıպ�������䣬��Ƭ���ڲ����ڵıպ����������
% 5.�޳��쳣ԭʼ����
% ===
% 6.��С���Բ���������Բ
% 7.��С��Ӿ�

rootdir = 'G:\workspace\���ݴ洢\��˿������ʵ��\Ƭ�̼��-2020.07.21';
% rootdir = 'G:\workspace\���ݴ洢\��˿������ʵ��\����ͼƬ����';
subdir = dir(rootdir);

% ����С��ֱ��趨��ͬ��ֵ
big_threshold = 800;
middle_threshold = 500;
small_threshold = 200;
detritus_threshold = 100;

% ���
% K = 0.0175;
% B = 0;

% ���
K = 0.0241;
B = 0;

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
    Atotal = cell(num_img,1); % cell�洢ÿ��ͼƬ��˿�����
    Ntotal = cell(num_img,1); % cell�洢ÿ��ͼƬ������
    
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
        
       %% ��ֵͼ��    
        gray = rgb2gray(I); % ͼ��ĻҶȴ���
        level = graythresh(gray);
        bw = imbinarize(gray,level);  % ͼ��Ķ�ֵ������
        bw1 = ~bw; % ���ÿ��Ƭ��
        bw2 =imfill(bw1,'holes');
        adjustImage = bwareaopen(bw2, threshold); % �޳�С���ͼ��
       
       %% ��ȡÿ����ͨ�������Ӧ����С���Բ
        regions = regionprops(adjustImage); % ������ͨ����region = [{area,centroid,boundingbox}]
        [boundaries,L] = bwboundaries(adjustImage); % ��ȡ�߽磬�˴����bw1���ն���䣬ʹ�߽���������ͨ����һ��
        
        [label,sheet_num] = bwlabel(adjustImage);
        Images = cell(1,sheet_num); % ÿ��Ƭ����������
        tarea = cell(1,sheet_num); % ÿ��Ƭ�����
        
        tcenter = cell(1,length(boundaries));
        tradius = zeros(1,length(boundaries));
        circlearea = zeros(1,length(boundaries));
        rrate = zeros(1,length(boundaries));
        
        tr = cell(1,sheet_num);
        for k = 1:sheet_num % ��ÿ��Ƭ�̷ֱ����
            tarea{k} = regions(k).Area;
            if tarea{k} > 1000000
                tarea{k} = []; %���Ƭ��84���趨Ƭ����ֵ������100��
                continue;
            end
            boundary = boundaries{k};
            [center,radius] = minboundcircle(boundary(:,2),boundary(:,1)); % ��С���Բ�뾶��Բ��
            tcenter{k} = center;
            tradius(k) = radius;
            circlearea(k) = pi * (radius.^2);
            
            rrate(k) = tarea{k} / circlearea(k);
            tr{k} = [tarea{k},rrate(k)];
        end
        Ntotal{j} = subdirfiles(j).name(1:end-4); 
        Atotal{j} = [tr{:}];
    end % j - ÿ��Ŀ¼�µ�ͼ��
    
    %% ͳ�ƽ����Ƭ����������Բ�����Բ����
    [curpath,~] = fileparts(mfilename('fullpath')); % ��ǰĿ¼
    dsave = strcat(curpath,'\Ƭ�����ͳ��.xlsx');
    writecell(Ntotal,dsave,'Sheet',i-2,'Range','A1'); % ��һ��д��ͼ���ļ���
    writecell(Atotal,dsave,'Sheet',i-2,'Range','B1'); % д��Ƭ�����
end %i - ��Ŀ¼