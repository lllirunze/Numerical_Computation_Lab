clc
clear
close

%% 初始化参数
a = -3;                     % 左边界
b = 3;                      % 右边界
c = 2;
d = 0.8;
e = -1;
f = 1.5;
n = 30;                     % [a,b]区间内平均划分数量
N = 100;
xi = linspace(a,b,n+1)';    % 在[a,b]之间共取n+1个数据点（中间有n-1个）
X = linspace(a,b,N+1)';
yi = c*sin(d*xi)+e*cos(f*xi);   % 准确的函数值
% yi(13) = yi(13)+0.2;        % 随机增加扰动（自己任意设定的，这句话怎么改都行）
plot(xi,yi,'bo');           % 绘制散点图
hold on
xlabel('x')
ylabel('y')
title('最小二乘法拟合f(x)=c*sin(d*x)+e*cos(f*x)')
hold on

%% 最小二乘法曲线拟合
p(:,1) = ones(n+1,1);
P(:,1) = ones(N+1,1);
for j=1:15
    p(:,2*j) = sin(j*xi);
    p(:,2*j+1) = cos(j*xi);
    P(:,2*j) = sin(j*X);
    P(:,2*j+1) = cos(j*X);
end
t = p\yi;
F = P*t;
plot(X,F,'g-');
hold on

%% 计算平均误差
m = 20;                     % 实验点的个数
testxi = round(linspace(2,N,m));   % 在[a,b]之间共取m个实验点
testyi = 1:20;
testlen = length(testxi);
testY1 = 1:20;
for i=1:testlen
    testyi(i) = c*sin(d*X(testxi(i)))+e*cos(f*X(testxi(i)));
    testY1(i) = F(testxi(i));
end
testerror = testyi-testY1;
% disp(testerror1);
Aver_error = 0;             % 平均误差
for i=1:testlen
    Aver_error = Aver_error+abs(testerror(i));
end
Aver_error = Aver_error/testlen;
disp(Aver_error);