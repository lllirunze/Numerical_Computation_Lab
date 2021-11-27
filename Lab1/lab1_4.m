clc
clear
close

%% 初始化参数
a = 0;                      % 左边界
b = 20;                     % 右边界
c = 2;
d = 0.8;
e = -1;
f = 1.5;
n = 10;                     % [a,b]区间内平均划分数量
xi = linspace(a,b,n+1);     % 在[a,b]之间共取n+1个数据点（中间有n-1个）
yi = c*sin(d*xi)+e*cos(f*xi);   % 准确的函数值
len = length(xi);           % 向量xi的长度
plot(xi,yi,'o');            % 绘制散点图
hold on
resX0 = xi(1):0.01:xi(end);
resY0 = c*sin(d*resX0)+e*cos(f*resX0);
h0 = plot(resX0,resY0,'k');
hold on
xlabel('x')
ylabel('y')
title('L(x)与f(x)的对比图像')
hold on
m = 20;                     % 实验点的个数
testxi = zeros(1,m);
testyi = zeros(1,m);
testY4 = zeros(1,m);
Aver_error = 0;             % 平均误差

%% 分段线性插值
for i=1:len-1
    resX4 = xi(i):0.01:xi(i+1);
    testxi(i*2-1) = xi(i);
    testxi(i*2) = (xi(i)+xi(i+1))/2;
    resY4 = yi(i)*(resX4-xi(i+1))/(xi(i)-xi(i+1))+yi(i+1)*(resX4-xi(i))/(xi(i+1)-xi(i));
    testyi(i*2-1) = c*sin(d*testxi(i*2-1))+e*cos(f*testxi(i*2-1));
    testY4(i*2-1) = yi(i)*(testxi(i*2-1)-xi(i+1))/(xi(i)-xi(i+1))+yi(i+1)*(testxi(i*2-1)-xi(i))/(xi(i+1)-xi(i));
    testyi(i*2) =  c*sin(d*testxi(i*2))+e*cos(f*testxi(i*2));
    testY4(i*2) = yi(i)*(testxi(i*2)-xi(i+1))/(xi(i)-xi(i+1))+yi(i+1)*(testxi(i*2)-xi(i))/(xi(i+1)-xi(i));
    h4 = plot(resX4,resY4,'r');
    hold on
end
legend([h0,h4],'f(x)=c*sin(d*x)+e*cos(f*x)','分段线性插值')

%% 计算平均误差
for i=1:m
   Aver_error = Aver_error+abs(testyi(i)-testY4(i)); 
end
Aver_error = Aver_error/m;
disp(Aver_error)