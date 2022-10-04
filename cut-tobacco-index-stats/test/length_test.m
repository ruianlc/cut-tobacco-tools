clc;clear;close all

dpath = 'C:\Users\Administrator\Documents\MATLAB\片烟烟丝指标\results\片烟模拟切丝结果\A.xlsx';
%dpath = 'C:\Users\Administrator\Documents\MATLAB\片烟烟丝指标\results\在线烟丝实验结果\#10.xlsx';

[~,Sheets,~] = xlsfinfo(dpath);
maxLens = cell(length(Sheets),1);
for ss = 1:length(Sheets)
    lenData = xlsread(dpath,ss);
    realLenData = lenData(:,4:end);
    
    [maxLenData] = max(max(realLenData));
    [rowId,colId] = find(realLenData == maxLenData);
    maxLens{ss} = [maxLenData,rowId,colId];
end

disp('bingo...')