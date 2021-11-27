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

%% 拉格朗日插值
u = [1 -xi(2)];             % u存放多项式系数
under = (xi(1) - xi(2));
for k=3:len
    v = [1 -xi(k)];         % v存放多次式系数
    u = conv(u,v);          % conv代表对两个多项式相乘
    under = under*(xi(1)-xi(k));
end
u = u/under;
u = u*yi(1);
temp1 = u;                  % 第一个项的结果

under = 1;
for k=2:len
    u = [1 -xi(1)];
    under = (xi(k)-xi(1));
    for j=2:len
        if j~=k
            v = [1 -xi(j)];
            u = conv(u,v);
            under = under*(xi(k)-xi(j));
        end
    end
    u = u/under;
    u = u*yi(k);
    temp2 = u;
    temp1 = temp1+temp2;
end
resX2 = xi(1):0.01:xi(end);
resY2 = polyval(temp1,resX2);
h2 = plot(resX2,resY2,'c');
hold on
legend([h0,h2],'f(x)=c*sin(d*x)+e*cos(f*x)','拉格朗日插值')

%% 计算平均误差
m = 20;                     % 实验点的个数
testxi = linspace(a,b,m);   % 在[a,b]之间共取m个实验点
testyi = c*sin(d*testxi)+e*cos(e*testxi);
testlen = length(testxi);
testY2 = polyval(temp1,testxi);
testerror2 = testyi-testY2;
% disp(testerror2);
Aver_error = 0;             % 平均误差
for i=1:testlen
    Aver_error = Aver_error+abs(testerror2(i));
end
Aver_error = Aver_error/testlen;
disp(Aver_error);