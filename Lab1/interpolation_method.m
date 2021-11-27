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
resX3 = xi(1):0.01:xi(end);
resY3 = subs(resY3,z,resX3);
h3 = plot(resX3,resY3,'g');
hold on

%% 分段线性插值
for i=1:len-1
    resX4 = xi(i):0.01:xi(i+1);
    resY4 = yi(i)*(resX4-xi(i+1))/(xi(i)-xi(i+1))+yi(i+1)*(resX4-xi(i))/(xi(i+1)-xi(i));
    h4 = plot(resX4,resY4,'r');
    hold on
end
    
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
    s1 = yi(1,i).*((resX5-xi(1,i+1)).^2).*(delta_h(1,i)+2.*(resX5-xi(1,i)))./(delta_h(1,i).^3);
    s2 = yi(1,i+1).*((resX5-xi(1,i)).^2).*(delta_h(1,i)+2.*(xi(1,i+1)-resX5))./(delta_h(1,i).^3);
    s3 = matrix2(1,i).*((resX5-xi(1,i+1)).^2).*(resX5-xi(1,i))./(delta_h(1,i).^2);
    s4 = matrix2(1,i+1).*((resX5-xi(1,i)).^2).*(resX5-xi(1,i+1))./(delta_h(1,i).^2);
    resY5 = s1+s2+s3+s4;
    h5 = plot(resX5,resY5,'b');
    hold on
end
    
%% 绘图
% 对不同曲线进行标注
% legend(h3,'牛顿插值')
% legend([h1,h2,h4,h5],'范德蒙德多项式插值','拉格朗日插值','分段线性插值','分段三次Hermite插值')
% legend([h0,h1,h2],'f(x)=c*sin(d*x)+e*cos(f*x)','范德蒙德多项式插值','拉格朗日插值')
legend([h0,h1,h2,h3,h4,h5],'f(x)=c*sin(d*x)+e*cos(f*x)','范德蒙德多项式插值','拉格朗日插值','牛顿插值','分段线性插值','分段三次Hermite插值')
hold on
