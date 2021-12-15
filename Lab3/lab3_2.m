clc
clear
close

a = 0;
b = 1;
n = 100;                % 将[a,b]划分为n等份
h = (b-a)/n;
syms x
y = (sin(x))/x;
temp = int(y,x,a,b);
I = sinint(1);          % I为sinx/x在[a,b]的定积分，即temp的表示
Sn = 0;
for k=0:(n-1)
    xk = a+k*h;
    xkhalf = xk+(h/2);
    xkplus = a+(k+1)*h;
    fx = integrand(xk);
    fxhalf = integrand(xkhalf);
    fxplus = integrand(xkplus);
    Sn = Sn+fx+4*fxhalf+fxplus;
end
Sn = Sn*h/6;            % 得到复合辛普森公式的值Sn
Rn = I-Sn;              % 得到余项Rn

function Y = integrand(X)   % 被积函数
    if X==0
        Y = 1;
    else
        Y = (sin(X))/X;
    end
end