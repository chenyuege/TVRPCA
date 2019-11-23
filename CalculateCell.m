function [MeanValue,StdValue] = CalculateCell(CellName)
load(CellName)
datalength  = length(MMPSNR);
methodwidth = length(MMPSNR{1});
MPSNR = zeros(datalength,methodwidth);
MSSIM = MPSNR;
ERGAS = MPSNR;
TIME = MPSNR;
MeanValue = zeros(4,methodwidth);
StdValue  = zeros(4,methodwidth);
for i =1:datalength
    MPSNR(i,:)= MMPSNR{i};
    MSSIM(i,:)= MMSSIM{i};
    ERGAS(i,:)= MMERGAS{i};
    TIME(i,:) = MMTIME{i};
end

MeanValue(1,:)= mean(MPSNR,1);
MeanValue(2,:)= mean(MSSIM,1);
MeanValue(3,:)= mean(ERGAS,1);
MeanValue(4,:)= mean(TIME,1);
StdValue(1,:) = std(MPSNR,1);
StdValue(2,:) = std(MSSIM,1);
StdValue(3,:) = std(ERGAS,1);
StdValue(4,:) = std(TIME,1);
% 取定四位小数
MeanValue = round(MeanValue,4);
StdValue = round(StdValue,4);
end



