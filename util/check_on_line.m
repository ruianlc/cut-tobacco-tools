function [lineFlag,residualsSum] = check_on_line(x,y,resiThreshold)

coefficients = polyfit(x, y, 1);
yFit = polyval(coefficients, x);
residualsSum = sum(abs(yFit - y));
if residualsSum < resiThreshold
        lineFlag = 1;
else
        lineFlag = 0;
end
