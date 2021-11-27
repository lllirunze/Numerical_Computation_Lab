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

%% 范德蒙德多项式插值
A = vander(xi);             % 利用范德蒙德矩阵生成插值多项式系数矩阵
resY1 = yi';
B = [A,resY1];              % 生成增广矩阵
C = rref(B);                % 化简，解方程组
D = C(:,end);
resX1 = xi(1):0.01:xi(end); % 绘制拟合曲线
resY1 = polyval(D,resX1);   % 写出拟合曲线方程
h1 = plot(resX1,resY1,'m');
hold on
legend([h0,h1],'f(x)=c*sin(d*x)+e*cos(f*x)','范德蒙德多项式插值')

%% 计算平均误差
m = 20;                     % 实验点的个数
testxi = linspace(a,b,m);   % 在[a,b]之间共取m个实验点
testyi = c*sin(d*testxi)+e*cos(e*testxi);
testlen = length(testxi);
testY1 = polyval(D,testxi);
testerror1 = testyi-testY1;
% disp(testerror1);
Aver_error = 0;             % 平均误差
for i=1:testlen
    Aver_error = Aver_error+abs(testerror1(i));
end
Aver_error = Aver_error/testlen;
disp(Aver_error);