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

rootdir = 'D:\workspace\���ݴ洢\���̾����˿\Bģ����˿ʵ��\������˿ʵ��\A';
%rootdir = 'D:\workspace\���ݴ洢\���̾����˿\Bģ����˿ʵ��\������˿ʵ��\B';
%rootdir = 'D:\workspace\���ݴ洢\���̾����˿\Bģ����˿ʵ��\������˿ʵ��\BС��';

subdir = dir(rootdir);

% ����С��ֱ��趨��ͬ��ֵ
big_threshold = 500;
middle_threshold = 200;
small_threshold = 100;
detritus_threshold = 30;

% ���
% K = 0.1051;
% B = 37.761;

% ����
K = 0.1314;
B = 0;

% 6���أ�0.8mm��8���أ�1.0mm��10���أ�1.2mm
PL = 8;%7.6104; % ������ʵ�ʳ��ȣ���˿���1mm����Ӧ����Ϊ7.6104����8��������Ϊһ����˿�Ŀ��
% Images = cell(1,sheet_num); % ÿ��ͼ������Ƭ����������

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
        source = imread(imgfile_path);
        I = source;
        rgb = reshape(I,size(I,1)*size(I,2),3);
        rgb(rgb(:,3)>90,:) = 0; % ����Bֵ����ȥ������
        rgb(rgb(:,3)~= 0,:) = 1;
        img = reshape(rgb,size(I,1),size(I,2),3);
        bw = img(:,:,1);
        BW = logical(bw);  
        adjustImage = bwareaopen(BW, threshold); % �޳�С���ͼ��

%         nhood = true(9);  % �˲�����
%         S = stdfilt(I,nhood);  % ��׼���˲�
%         gray = mat2gray(S);  % �ҶȻ�
%         umb = rgb2gray(gray); % ͼ��ĻҶȴ���
%         level = graythresh(umb);
%         bw = imbinarize(umb,level);  % ͼ��Ķ�ֵ������
%         bw1 = ~bw; % ���ÿ��Ƭ��
%         bw2 = imfill(bw1,'holes');
%         adjustImage = bwareaopen(bw2, threshold); % �޳�С���ͼ��
%         figure
%         imshow(adjustImage)
        
        %% ��ȡ����߽�
        [boundaries,L] = bwboundaries(adjustImage);
        area = zeros(1,length(boundaries));  % ����Ƭ�����
        tarea = zeros(1,length(boundaries)); % ����Ƭ�����Բ���
        rrate = zeros(1,length(boundaries)); % ����Ƭ��Բ����
   
        %% ��ȡÿ����ͨ����
        regions = regionprops(adjustImage); % ������ͨ����region = [{area,centroid,boundingbox}]
        [label,sheet_num] = bwlabel(adjustImage);

        Limage = cell(sheet_num,1);% ÿ��ͼ��Ƭ��ָ�꣺�����Բ���ʡ���˿����
        
        for s = 1:sheet_num % ��ÿ��Ƭ�̷ֱ����
            bw = label;
            pos = regions(s).BoundingBox;
            
            r1 = round(pos(2));
            c1 = round(pos(1));
            w = pos(3);
            h = pos(4);
            r2 = r1+h-1;
            c2 = c1+w-1;
            bw(bw ~= s) = 0;  % ����������0
            bw(bw == s) = 1;  % ������i��������1
            subImg = bw(r1:r2,c1:c2);  % ÿ��Ƭ�̵���������
            area(s) = bwarea(subImg);  % ����Ƭ�����
            
%             figure
%             imshow(Images{s})
            
            boundary = boundaries{s};
            [center,radius] = minboundcircle(boundary(:,2),boundary(:,1)); % ��С���Բ�뾶��Բ��
            
            % �������Բ
            % viscircles(center,radius);
            % tcenter{i} = center;
            % tradius(i) = radius;
            % tarea(s) = pi * (radius.^2);
            % rrate(s) = area(s) / tarea(s);
            
            %% Ƭ���и�
            % �Ե�����ͨ�������ģ����˿
            % �Կ��wΪ�ȷֶ�����10��������Ϊ3mm����Ϊ1����˿�Ŀ�ȣ�����ÿ����ͨ�����и������˿����
            stick_num = floor(w/PL);
            Lsheet = cell(1,stick_num); % һ��Ƭ�̵�������˿����

            %% Ƭ���и���˿����
            for t = 1:stick_num

                subStick = subImg(:,(t-1)*PL+1:t*PL);
                
                % �жϴ˴�����ͨ����������и�������˿
                stick_regions = regionprops(subStick);
                stick_area = bwlabel(subStick);
                num_stick_regions = length(stick_regions);
                sStick = zeros(1,num_stick_regions); % �и���˿�ĳ��ȣ����ڶ�˿��
                
                if num_stick_regions == 1 % δ������˿
                    skL = sketelon(subStick); % ��ȡ�Ǽ�
                    [ROWS,COLUMNS] = find(skL == 1);  % ��ȡ�Ǽܵ�����
                    
                    if size(ROWS,1) <= 20 % �޳�������С���
                        continue;
                    end
                    DDist = count_len(ROWS,COLUMNS);
                    DDist(DDist > 1.5) = []; % �޳����Ե�غ���˿���쳣ֵ
                    Pixel = sum(DDist);
                    Lsheet{t} = K * Pixel + B;
                else % ������˿
                    
                    for p = 1:num_stick_regions % ÿ����˿�����Ķ�˿����
                        bw_stick = stick_area;
                        pos_stick = stick_regions(t).BoundingBox; % ��Ӿ��ε�x��y�߽磨�½�����䣩��pos = [x,y,width,height]
                        
                        r1 = round(pos_stick(2));
                        c1 = round(pos_stick(1));
                        w = pos_stick(3);
                        h = pos_stick(4);
                        r2 = r1+h-1;
                        c2 = c1+w-1;
                        bw_stick(bw_stick ~= p) = 0;
                        bw_stick(bw_stick == p) = 1;
                        subsubStick = bw_stick(r1:r2,c1:c2);
                        
                        sskL = sketelon(subsubStick); % ��ȡ�Ǽ�
                        [ROWS,COLUMNS] = find(sskL == 1);  % ��ȡ�Ǽܵ�����
                        
                        if size(ROWS,1) <= 20 % �޳�������С���
                            continue;
                        end
                        DDist = count_len(ROWS,COLUMNS);
                        DDist(DDist > 1.5) = []; % �޳����Ե�غ���˿���쳣ֵ
                        Pixel = sum(DDist);
                        sStick(p) = K * Pixel + B;
                    end
                    Lsheet{t} = sStick(:);
                end
            end
            Limage{s} = [area(s),Lsheet{:}];
        end
        Ltotal = [Ltotal;Limage];
    end % j - ÿ��Ŀ¼�µ�ͼ��
    
    %% ͳ�ƽ����Ƭ����������Բ�����Բ����
    [curpath,~] = fileparts(mfilename('fullpath')); % ��ǰĿ¼
    dsave = strcat(curpath,'\results\Ƭ��ģ����˿���\','A.xlsx');
    
    stitle = {'���','Ƭ����������أ�','��˿���ȣ�mm��'};
    xuhao = 1:length(Ltotal);
    writecell(stitle,dsave,'Sheet',i-2,'Range','A1'); % д���ͷ
    writematrix(xuhao',dsave,'Sheet',i-2,'Range','A2'); % д�����
    writecell(Ltotal,dsave,'Sheet',i-2,'Range','B2'); % д��Ƭ�����
end %i - ��Ŀ¼