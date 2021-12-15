clc
clear
close

a = 0;
b = 1;
syms x
y = (sin(x))/x;
temp = int(y,x,a,b);
real = sinint(1);       % real为sinx/x在[a,b]的定积分，即temp的表示
h = b-a;                % 最大步长
n = 10;                 % 龙贝格算法的精度
epsilon = 0.000001;     % 预先给定的精度
T = zeros(n,n);
fa = integrand(a);
fb = integrand(b);
T(0+1,0+1) = (fa+fb)*h/2;
k = 0;                  % T表的行
m = 0;                  % T表的列
step = 1;               % 步数
err1 = 1;               % 初始差值
while (err1>epsilon)&&(k<n)
    k = k+1;
    h = h/2;
    sum = 0;
    for  p=1:step
        xhalf = a+h*(2*p-1);
        fhalf = integrand(xhalf);
        sum = sum+fhalf;
    end
    T(k+1,1) = T(k,1)/2+h*sum;  % 公式(4,1)计算每行的第一列值
    step = step*2;      % 更新步数
    for m=1:k
        % 公式(4,11)计算每一行的值
        T(k+1,m+1) = T(k+1,m)+(T(k+1,m)-T(k,m))/(4^m-1);
    end
    err1 = abs(T(k,k)-T(k+1,k+1));
    % 计算误差error
    error = abs(T(k+1,k+1)-real);
end
I = T(k+1,k+1);         % 终止循环，最后的值约等于I

function Y = integrand(X)   % 被积函数
    if X==0
        Y = 1;
    else
        Y = (sin(X))/X;
    end
end