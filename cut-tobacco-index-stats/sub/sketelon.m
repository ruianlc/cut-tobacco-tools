function skL = sketelon(bws)
%% �Ǽ���ȡ����
%% input:
% bws����������Ƭ�̶�ֵͼ�񣨾����и��
%% output:
% skL��������ȡ�ĹǼ�
%% ����ʼ
sk = bwmorph(bws,'skel',Inf);  % ��ȡ�Ǽ�
bm = ones(3); % ����3*3��ȫΪ1�ľ���
bm(2,2) = 0;  % ���м����ֵ����Ϊ0
[ro,la] = size(bws); % ��ֵͼ��Ĵ�С
ze1 = zeros(ro,1);  % ����ͬ��ֵͼ��ͬ������ȫ0������
ze2 = zeros(1,la+2); % ����һ�ж�ֵͼ����+2��������
sk1 = [ze1,sk,ze1];  % ���Ǽ�ͼ���������߸���һ��0����
sk2 = [ze2;sk1;ze2];  % ���Ǽ����¸���һ��0����
sumsk = zeros(ro,la);  % ����ȫ0��ͬ��ֵͼ���С��0����
%% Ѱ�ҽڵ�
for i = 2:ro
    for j = 2:la
        SK = sk2(i-1:i+1,j-1:j+1);  % ÿ����ȡ3*3�ľ���
        SK1 = SK .* bm;  % ���þ���������ֵ��0
        su = sum(sum(SK1));  % ��þ���ĺ�
        if (su > 2) && (sk2(i,j) ~= 0)  % ����ʹ���2���Ҹõ�ֵ������0
            sumsk(i-1,j-1) = su;  % �ڵ��ֵ��ͼ
        end
    end
end
[front1,fn1] = bwlabel(sumsk); % ͼ����
skk = sk;  % �Ǽ�ͼ��
szo = zeros(ro,la);  % ȫ���ֵͼ��
for m = 1:fn1
    sk0 = skk;  % �Ǽ�ͼ
    sk00 = skk;  % �Ǽ�ͼ
    sk0(find(front1==m)) = 0;  % ���Ǽ�ͼ���иýڵ�ɾ��
    [front2,fn2] = bwlabel(sk0); % ���ɾ���ýڵ���ͼ��
    if fn2 < 3  % �����֧С��3������������
        continue
    end
    farea = front2(:);  % ������ת������
    [fcounts,~] = histc(farea,1:fn2);  % ͳ�Ƹ���֧�ĳ���
    [~,tn] = sort(fcounts); % �Գ��Ƚ���������С����
    %[~,st]=min(fcounts);
    for i = 1:length(tn)-2
        sk00(find(front2 == tn(i))) = 0;  % �������������֧
    end
    ss = skk - sk00;  % ��֧ͼ��
    szo = ss + szo;  % ÿ���ڵ��֧
end
szo(find(szo ~= 0)) = 1;
skL = skk - szo;  % ȥ����֧��Ľڵ�ͼ��
skL = logical(skL);
skL = bwmorph(skL,'spur',1);  % ���޳�С�ķ�֧
