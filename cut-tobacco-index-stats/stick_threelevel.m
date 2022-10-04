clc;clear;close all

rootdir = 'G:\workspace\���ݴ洢\������˿���\����ʵ��-��˿���ȼ��';

subdir = dir(rootdir);

% ���ж���ֱ��趨��ͬ��ֵ
long_threshold = 300;
middle_threshold = 200;
short_threshold = 60;
detritus_threshold = 10;

for i = 1:length(subdir)
    
    % ����'.'��'..'Ŀ¼
    if(strcmp(subdir(i).name,'.')||...
            strcmp(subdir(i).name,'..')||...
            ~subdir(i).isdir())
        continue;
    end
    
    subsubdir =dir(strcat(rootdir,'\',subdir(i).name));
    for j = 1:length(subsubdir)
        
        % ����'.'��'..'Ŀ¼
        if(strcmp(subsubdir(j).name,'.')||...
                strcmp(subsubdir(j).name,'..')||...
                ~subsubdir(j).isdir())
            continue;
        end
        
        % ���ļ������Һ�׺Ϊbmp���ļ�
        subsubdirpath = fullfile( rootdir,subdir(i).name, subsubdir(j).name);
        subsubdirfiles = dir(strcat(subsubdirpath,'\','*.bmp'));
        num_img = length(subsubdirfiles); % ÿ�����ļ���ͼ������
        Ltotal = cell(num_img,1); % cell�洢ÿ��ͼƬ��˿�ĳ���
        Ntotal = cell(num_img,1); % cell�洢ÿ��ͼƬ������
        if(isequal(num_img,0))
            return;
        end
        
        if contains(subsubdir(j).name,'��')
            threshold = long_threshold;
        elseif contains(subsubdir(j).name,'��')
            threshold = middle_threshold;
        elseif contains(subsubdir(j).name,'��')
            threshold = short_threshold;
        elseif contains(subsubdir(j).name,'��')
            threshold = detritus_threshold;
        else
            threshold = middle_threshold;
        end
        
        for s = 1:num_img
            
            imgfile_path = fullfile(subsubdirpath, subsubdirfiles(s).name);
            I = imread(imgfile_path);
%             figure
%             imshow(I)
%             titleString = strcat('��',num2str(i),'������ -- ','�ļ��У�',subsubdir(j).name,' -- ��',num2str(s),'��ͼ��');
%             title(titleString)
            
            source = I(1:end-100,200:end-200,:); %��ȥ���½���ӰӰ��
%             figure
%             imshow(source)
%             title('����ͼ��')

            rgb = reshape(source,size(source,1)*size(source,2),3);
            rgb(rgb(:,3)>240,:) = 0; % ȥ��������Bֵ����240
            rgb(rgb(:,3)~= 0,:) = 1;
            img = reshape(rgb,size(source,1),size(source,2),3);
            bw0 = img(:,:,1);
            bw = logical(bw0);
%             figure
%             imshow(bw)
%             title('��ֵ��ͼ��')
            
            BW2 = imfill(bw,'holes');  % ���ͼ��
%             figure
%             imshow(BW2)
%             title('��ֵ�����ͼ��')
            
            adjustImage = bwareaopen(BW2, threshold); % �޳�С���ͼ��
%             figure
%             imshow(adjustImage)
%             title('��ֵ���޳�С���ͼ��')
%             hold on
            
            regions = regionprops(adjustImage); % ������ͨ����
            [labels,Num] = bwlabel(adjustImage);
            
            % ��ȡ������˿ͼ��
%             K = 0.3262;
%             B = 0.7791;
            K = 0.1549; % �������غ���ת��
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
        
                skL = sketelon(stick); % ��ȡ�Ǽ�
%                 figure
%                 imshow(skL)
                
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
            
            % ������
            Ntotal{s} = strcat(subsubdir(j).name,subsubdirfiles(s).name(end-6:end)); % ÿ��ͼ�����ƣ�ȥ���ļ���׺��
            Ltotal{s} = Lens;  % ÿ��ͼ���������˿����
        end
        [curpath,~] = fileparts(mfilename('fullpath'));
        desPath = strcat(curpath,'\reults\',subdir(i).name,'.xlsx');
        writecell(Ntotal,desPath,'Sheet',j-2,'Range','A1'); % д����˿����
        writecell(Ltotal,desPath,'Sheet',j-2,'Range','B1'); % д����˿����
    end
end