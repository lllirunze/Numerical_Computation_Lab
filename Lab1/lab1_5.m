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

%% 分段三次Hermite插值
syms t
defivative_a = limit(c*sin(d*t)+e*cos(f*t),a);  % 函数f(x)在a点的一阶导数 
defivative_b = limit(c*sin(d*t)+e*cos(f*t),b);  % 函数f(x)在b点的一阶导数
[~,number] = size(xi);                          % 获取输入xi的大小1*number
delta_h = zeros(1,number-1);                    % 给delta_h分配数组大小1*(number-1)，并全部初始化为0
delta_f = zeros(1,number-1);                    % 给delta_f分配数组大小1*(number-1)，并全部初始化为0
lambda_ = zeros(1,number-2);                    % 给lambda_分配数组大小1*(number-2)，并全部初始化为0
miu_ = zeros(1,number-2);                       % 给miu_分配数组大小1*(number-2)，并全部初始化为0
sigma_ = zeros(1,number-2);                     % 给sigma_分配数组大小1*(number-2)，并全部初始化为0
X = zeros(number-2,number-2);                   % 初始化系数矩阵X，(n-2)*(n-2)
Y = zeros(number-2,1);                          % 初始化系数矩阵Y，(n-2)*1
% 计算delta_h, delta_f的值
for i=1:(number-1)
    delta_h(i) = xi(i+1)-xi(i);
    delta_f(i) = (yi(i+1)-yi(i))/delta_h(i);
end
% 计算lambda_, miu_, sigma_的值
for i=1:(number-2)
    lambda_(1,i) = delta_h(1,i+1)/(delta_h(1,i+1)+delta_h(1,i));
    miu_(1,i) = 1-lambda_(1,i);
    sigma_(1,i) = 3*(lambda_(1,i)*delta_f(1,i)+miu_(1,i)*delta_f(1,i+1));
end
% 分三种情况，i=1, i=2:n-2, i=n-1
X(1,1) = 2;
X(1,2) = miu_(1,1);
Y(1,1) = sigma_(1,1)-lambda_(1,1)*defivative_a;
for i=2:number-3
    Y(i,1) = sigma_(1,i);
    X(i,i-1) = lambda_(1,i);
    X(i,i) = 2;
    X(i,i+1) = miu_(1,i);
end
X(number-2,number-3) = lambda_(1,number-2);
X(number-2,number-2) = 2;
Y(number-2,1) = sigma_(1,number-2)-miu_(1,number-2)*defivative_b;
% 计算A的逆，A*B
matrix1 = X\Y;
matrix2 = zeros(1,number);
matrix2(1,1) = defivative_a;
matrix2(1,number) = defivative_b;
for i=2:number-1
    matrix2(1,i) = matrix1(i-1,1);
end
for i=1:number-1
    % 获取相邻两点之间的插值函数
    resX5 = linspace(xi(1,i),xi(1,i+1));
    testxi(i*2-1) = xi(i);
    testxi(i*2) = (xi(i)+xi(i+1))/2;
    s1 = yi(1,i).*((resX5-xi(1,i+1)).^2).*(delta_h(1,i)+2.*(resX5-xi(1,i)))./(delta_h(1,i).^3);
    s2 = yi(1,i+1).*((resX5-xi(1,i)).^2).*(delta_h(1,i)+2.*(xi(1,i+1)-resX5))./(delta_h(1,i).^3);
    s3 = matrix2(1,i).*((resX5-xi(1,i+1)).^2).*(resX5-xi(1,i))./(delta_h(1,i).^2);
    s4 = matrix2(1,i+1).*((resX5-xi(1,i)).^2).*(resX5-xi(1,i+1))./(delta_h(1,i).^2);
    resY5 = s1+s2+s3+s4;
    testyi(i*2-1) = c*sin(d*testxi(i*2-1))+e*cos(f*testxi(i*2-1));
    s1_ = yi(1,i).*((testxi(i*2-1)-xi(1,i+1)).^2).*(delta_h(1,i)+2.*(testxi(i*2-1)-xi(1,i)))./(delta_h(1,i).^3);
    s2_ = yi(1,i+1).*((testxi(i*2-1)-xi(1,i)).^2).*(delta_h(1,i)+2.*(xi(1,i+1)-testxi(i*2-1)))./(delta_h(1,i).^3);
    s3_ = matrix2(1,i).*((testxi(i*2-1)-xi(1,i+1)).^2).*(testxi(i*2-1)-xi(1,i))./(delta_h(1,i).^2);
    s4_ = matrix2(1,i+1).*((testxi(i*2-1)-xi(1,i)).^2).*(testxi(i*2-1)-xi(1,i+1))./(delta_h(1,i).^2);
    testY4(i*2-1) = s1_+s2_+s3_+s4_;
    testyi(i*2) =  c*sin(d*testxi(i*2))+e*cos(f*testxi(i*2));
    s1__ = yi(1,i).*((testxi(i*2)-xi(1,i+1)).^2).*(delta_h(1,i)+2.*(testxi(i*2)-xi(1,i)))./(delta_h(1,i).^3);
    s2__ = yi(1,i+1).*((testxi(i*2)-xi(1,i)).^2).*(delta_h(1,i)+2.*(xi(1,i+1)-testxi(i*2)))./(delta_h(1,i).^3);
    s3__ = matrix2(1,i).*((testxi(i*2)-xi(1,i+1)).^2).*(testxi(i*2)-xi(1,i))./(delta_h(1,i).^2);
    s4__ = matrix2(1,i+1).*((testxi(i*2)-xi(1,i)).^2).*(testxi(i*2)-xi(1,i+1))./(delta_h(1,i).^2);
    testY4(i*2) = s1__+s2__+s3__+s4__;
    h5 = plot(resX5,resY5,'b');
    hold on
end
legend([h0,h5],'f(x)=c*sin(d*x)+e*cos(f*x)','分段三次Hermite插值')

%% 计算平均误差
for i=1:m
   Aver_error = Aver_error+abs(testyi(i)-testY4(i)); 
end
Aver_error = Aver_error/m;
disp(Aver_error)