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

% rootdir = 'G:\workspace\���ݴ洢\������˿���\��˿������ʵ��\Ƭ�̼��-2020.07.21';
% rootdir = 'G:\workspace\���ݴ洢\������˿���\��˿������ʵ��\����ͼƬ����';
rootdir = 'G:\workspace\���ݴ洢\������˿���\������˿ʵ��\A';
subdir = dir(rootdir);

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
    Ltotal = cell(num_img,1); % cell�洢ÿ��ͼƬ��˿�ĳ���
    Ntotal = cell(num_img,1); % cell�洢ÿ��ͼƬ������
    
    if(isequal(num_img,0))
        return;
    end
    
    for j = 1:num_img
        
        imgfile_path = fullfile(subdirpath, subdirfiles(j).name);
        I = imread(imgfile_path);
        
       %% ��ֵͼ��    
        gray = rgb2gray(I); % ͼ��ĻҶȴ���
        level = graythresh(gray);
        bw = imbinarize(gray,level);  % ͼ��Ķ�ֵ������
        bw1 = ~bw; % ���ÿ��Ƭ��
        adjustImage = bwareaopen(bw1, 300); % �޳�С���ͼ��
       
        % ��ȡÿ����ͨ����
        regions = regionprops(adjustImage); % ������ͨ����region = [{area,centroid,boundingbox}]
        [area,sheet_num] = bwlabel(adjustImage);
        tpos=regionprops(area,'BoundingBox');
        tcentroid = regionprops(area,'Centroid');
        
        % �����ע��ͨ����
        % for i = 1:sheet_num % ÿ��Ƭ��
        %     rectangle('position',tpos(i).BoundingBox,'edgecolor','r', 'LineWidth',1.5);
        %     line([tcentroid(i,1).Centroid(1,1), tcentroid(i,1).Centroid(1,1)],[tcentroid(i,1).Centroid(1,2)-35, tcentroid(i,1).Centroid(1,2)+35],'Color','b', 'LineWidth',2);
        %     text(tcentroid(i,1).Centroid(1,1),tcentroid(i,1).Centroid(1,2), num2str(i),'Color', 'r');
        % end
        
        K = 0.1314;
        B = 0;
        
        PL = 7.6104; % ������ʵ�ʳ��ȣ���˿���1mm����Ӧ����Ϊ7.6104
        Images = cell(1,sheet_num);
        LL = cell(1,sheet_num);% ÿ��ͼ��������˿����
       
        % ���ÿ����ͨ�������Ӿأ��õ�����
        for k = 1:sheet_num % ��ÿ��Ƭ�̷ֲ�����
            bw = area;
            pos = regions(k).BoundingBox; % ��Ӿ��ε�x��y�߽磨�½�����䣩��pos = [x,y,width,height]
            
            r1 = round(pos(2));
            c1 = round(pos(1));
            w = pos(3);
            h = pos(4);
            r2 = r1+h-1;
            c2 = c1+w-1;
            bw(bw ~= k) = 0;  % ����������
            bw(bw == k) = 1;  % ������k��������1
            Images{k} = bw(r1:r2,c1:c2);  % ÿ��Ƭ�̵���������
            
            % �Ե�����ͨ�������ģ����˿���Կ��wΪ�ȷֶ�����10��������Ϊ3mm����Ϊ1����˿�Ŀ�ȣ�����ÿ����ͨ�����и������˿����
            stick_num = floor(w/PL);
            Stick = cell(1,stick_num); % ÿ����˿
            skL = cell(1,stick_num); % ÿ����˿�Ǽ�
            L = cell(1,stick_num); % һ��Ƭ�̵�������˿����
            
            for t = 1:stick_num
                
                % ��ȡÿ����˿��Images{k}-->Stick{t}
                Stick{t} = Images{k}(:,(t-1)*PL+1:t*PL);
                
              %% TODO �жϴ˴�����ͨ����������и�������˿
                stick_regions = regionprops(Stick{t}); 
                stick_area = bwlabel(Stick{t});
                num_stick_regions = length(stick_regions);
                sStick = zeros(1,num_stick_regions);
                
                if num_stick_regions == 1 % δ������˿
                    skL{t} = sketelon(Stick{t}); % ��ȡ�Ǽ�
                    [ROWS,COLUMNS] = find(skL{t} == 1);  % ��ȡ�Ǽܵ�����
                    
                    if size(ROWS,1) <= 20 % �޳�������С���
                        continue;
                    end
                    DDist = count_len(ROWS,COLUMNS);
                    Pixel = sum(DDist);
                    L{t} = K * Pixel + B;
                else % ������˿
                    
                    for p = 1:num_stick_regions % ÿ����˿�����Ķ�˿����
                        bw_stick = stick_area;
                        pos_stick = stick_regions(k).BoundingBox; % ��Ӿ��ε�x��y�߽磨�½�����䣩��pos = [x,y,width,height]
                        
                        r1 = round(pos_stick(2));
                        c1 = round(pos_stick(1));
                        w = pos_stick(3);
                        h = pos_stick(4);
                        r2 = r1+h-1;
                        c2 = c1+w-1;
                        bw_stick(bw_stick ~= p) = 0;  
                        bw_stick(bw_stick == p) = 1; 
                        subStick{p} = bw_stick(r1:r2,c1:c2);  
                        
                        sskL{p} = sketelon(subStick{p}); % ��ȡ�Ǽ�
                        [ROWS,COLUMNS] = find(sskL{p} == 1);  % ��ȡ�Ǽܵ�����
                        
                        if size(ROWS,1) <= 20 % �޳�������С���
                            continue;
                        end
                        DDist = count_len(ROWS,COLUMNS);
                        Pixel = sum(DDist);
                        sStick(p) = K * Pixel + B;
                    end
                    L{t} = sStick(1,:);
                end
            end % t - ÿ��Ƭ�̵���˿
            LL{k} = [L{:}];
        end %k - ÿ��ͼ���е�Ƭ�� 
        Ntotal{j} = subdirfiles(j).name(1:end-4); 
        Ltotal{j} = [LL{:}];
    end % j - ÿ��Ŀ¼�µ�ͼ��
    
    [curpath,~] = fileparts(mfilename('fullpath')); % ��ǰĿ¼
    dsave = strcat(curpath,'\Ƭ���и���˿����ͳ��.xlsx');
    writecell(Ntotal,dsave,'Sheet',i-2,'Range','A1'); % ��һ��д��ͼ���ļ���
    writecell(Ltotal,dsave,'Sheet',i-2,'Range','B1'); % д����˿����
end %i - ��Ŀ¼