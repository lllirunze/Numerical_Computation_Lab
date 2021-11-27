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

%% 牛顿插值
matrix0 = zeros(len,len);
% 对差商表第一列赋值
for i=1:len
    matrix0(i) = yi(i);
end
% 求差商表，从0阶开始，但从1维开始存储
for i=2:len
    for j=i:len
        matrix0(j,i) = (matrix0(j,i-1)-matrix0(j-1,i-1))/(xi(j)-xi(j-i+1));
    end
end
% disp('差商表：');
% disp(matrix0);
% 求牛顿插值多项式
z = sym('z');
resY3 = 0;
for k=2:n
    temp3 = 1;
    for j=1:k-1
        temp3 = temp3*(z-xi(j));
        % disp(temp3)
    end
    resY3 = matrix0(k,k)*temp3+resY3;
    % disp(resY3)
end
resY3 = matrix0(1,1)+resY3;
testY3 = resY3;
resX3 = xi(1):0.01:xi(end);
resY3 = subs(resY3,z,resX3);
h3 = plot(resX3,resY3,'g');
hold on
legend([h0,h3],'f(x)=c*sin(d*x)+e*cos(f*x)','牛顿插值')

%% 计算平均误差
m = 20;                     % 实验点的个数
testxi = linspace(a,b,m);   % 在[a,b]之间共取m个实验点
testyi = c*sin(d*testxi)+e*cos(e*testxi);
testlen = length(testxi);
testY3 = subs(testY3,z,testxi);
testerror3 = testyi-testY3;
% disp(testerror3);
Aver_error = 0;             % 平均误差
for i=1:testlen
    Aver_error = Aver_error+abs(testerror3(i));
end
Aver_error = Aver_error/testlen;
disp(Aver_error);