clc
clear
close

syms x
A = [   30,   33,  -43,  -11,  -38,  -29,   37,   28,   23;
      -480, -523,  664,  128,  621,  480, -618, -489,  329;
        60,  266,-1862,-1991,  464,  546, -968,-1567, 1652;
       540,  624, -782,  290, -893,  123,  567,    5, -122;
      -450, -675, 2245, 2326,-1512, 1230, -822,  129, -189;
      -300, -120,-1114,-1295, 1946,  302, -376,-1540, -609;
      1080,  998,  508, 2460,-1628,-1358, 2896, 2828,-2002;
     -1080,-1408, 3340, 2267,   21,-1202,  866,-2690,-1351;
      -300, -435, 1594, 1685,  340, 2279,  -27, 2917,-2336];
b = [  188;
     -3145;
     -4994;
       680;
      7845;
      1876;
      9712;
    -11599;
     10127];
% LU分解法
[L,U,x] = LU(A,b);

syms xrand
n = 20;
min = -1000; % 产生随机数的最小值
max = 1000;  % 产生随机数的最大值
Arand = round(min+(max-min)*rand(n,n));
brand = round(min+(max-min)*rand(n,1));
[Lrand,Urand,xrand] = LU(Arand,brand);

function [L,U,x] = LU(a,b)
    n = length(a);
    y = zeros(n,1);
    x = zeros(n,1);
    L = zeros(n,n);
    U = zeros(n,n);
    for i=1:n
        L(i,i) = 1;
    end
    for k=1:n
        for j=k:n
            U(k,j) = a(k,j)-sum(L(k,1:k-1)*U(1:k-1,j));
        end
        for i=k+1:n
            L(i,k) = (a(i,k)-sum(L(i,1:k-1)*U(1:k-1,k)))/U(k,k);
        end
    end
    y(1) = b(1)/L(1,1);
    for k=2:n
        for j=1:k-1
            y(k) = (b(k)-sum(L(k,1:j)*y(1:j)))/L(k,k);
        end
    end
    x(n) = y(n)/U(n,n);
    for i=n-1:-1:1
        x(i) = (y(i)-sum(U(i,i+1:n)*x(i+1:n)))/U(i,i);
    end
end
