clc;clear;close all
%%
% --- ˵�� ---
% ��ȡ����ļ����µ�Ҷ˿ͼ�񣬼���ÿ��ͼ���е�Ҷ˿����
% 
% Programmer: Robin An, 2021.07.26
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rootdir = 'G:\workspace\���ݴ洢\������˿���\��˿������ʵ��\Ҷ˿���-2021.07.22';
%rootdir = 'G:\workspace\���ݴ洢\��˿������ʵ��\����ͼƬ����';
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
        
        % �Ե���ͼ��Ĵ���
        source = I(80:end-200,80:end-150,:);  % ͼ���и�
        % ͼ���˲�
        nhood = true(9);  % �˲�����
        S = stdfilt(source,nhood);  % ��׼���˲�
        gray = mat2gray(S);  % �ҶȻ�
        umb = rgb2gray(gray); % ͼ��ĻҶȴ���
        bw = imbinarize(umb,0.15);  % ͼ��Ķ�ֵ������im2bw��rgbͼ��Ϊ���룬�ٷ��Ѳ��Ƽ�ʹ�ã�imbinarize�ԻҶ�ͼ��Ϊ����
        BW2 = imfill(bw,'holes');  % ���ͼ��
        adjustImage = bwareaopen(BW2, 100); % �޳�С���ͼ��
        contour = bwperim(adjustImage); % ��ȡ����
        regions = regionprops(adjustImage); % ������ͨ����
        area = bwlabel(adjustImage);
        % [area,Num] = bwlabel(adjustImage);
        
        
        % ��ȡÿ����˿ͼ��
        % K=0.1007;
        % B=1.641;
        K = 1;
        B = 0;
        Num = length(regions);
        L = zeros(1,Num); % ����ͼ����������˿�ĳ���
                
        for idx = 1:Num % idxΪ����ͼ����ÿ����˿�����
            bw = area;
            pos = regions(idx).BoundingBox; % ��С��Ӿ��ε�x��y�߽磨�½�����䣩��pos = [x,y,width,height]
            r1 = round(pos(2));
            c1 = round(pos(1));
            w = pos(3);
            h = pos(4);
            r2 = r1+h-1;
            c2 = c1+w-1;
            bw(bw ~= idx) = 0;  % ����������
            bw(bw == idx) = 1;  % ������i��������1
            Images{idx} = bw(r1:r2,c1:c2);  % ������˿ͼ��
            figure
            imshow(Images{idx})
            % Images{idx} = bw(r1:r2,c1:c2,:); 
            contour = bwperim(Images{idx});
            I1 = bwmorph(Images{idx},'skel',inf);
            [endpoint_index, nodepoint_index, rows, columns] = endnode(I1);
            [x1,y1] = find(contour == 1);
            [x2,y2] = find(I1 == 1);
            skL{idx} = sketelon(Images{idx}); % ��ȡ�Ǽ�
            [ROWS,COLUMNS] = find(skL{idx} == 1);  % ��ȡ�Ǽܵ�����
            DDist = count_len(ROWS,COLUMNS);
            Pixel = sum(DDist);
            L(idx) = K * Pixel + B - 3;
        end
        
        % ������
        Ntotal{j} = subdirfiles(j).name(1:end-4); % ÿ��ͼ�����ƣ�ȥ���ļ���׺��
        Ltotal{j} = L;  % ÿ��ͼ���������˿����
    end
    [curpath,~] = fileparts(mfilename('fullpath')); % ��ǰĿ¼
    dsave = strcat(curpath,'\��˿����ͳ��.xlsx');
    writecell(Ntotal,dsave,'Sheet',i-2,'Range','A1'); % ��һ��д��ͼ���ļ���
    writecell(Ltotal,dsave,'Sheet',i-2,'Range','B1'); % д����˿����
end