clc
clear
close

% 输入量：
% A: 线性方程组的系数矩阵（n*n，非奇异）
% b: 方程组右边的常数项列向量
% n: 方程组维数
% x0: 初始值
% N:  最大迭代次数
% 
% 迭代终止标准：
% 已达精度上限值或者到达最大迭代次数
% 
% 输出量：
% x：线性方程组的解

A = [ 31, 13,  0,  0,  0,-10,  0,  0,  0;
     -13, 35, -9,  0,-11,  0,  0,  0,  0;
       0, -9, 31,-10,  0,  0,  0,  0,  0;
       0,  0,-10, 79,-30,  0,  0,  0, -9;
       0,  0,  0,-30, 57, -7,  0, -5,  0;
       0,  0,  0,  0, -7, 47,-30,  0,  0;
       0,  0,  0,  0,  0,-30, 41,  0,  0;
       0,  0,  0,  0, -5,  0,  0, 27, -2;
       0,  0,  0, -9,  0,  0,  0, -2, 29];
b = [-15;
      27;
     -23;
       0;
     -20;
      12;
      -7;
       7;
      10];
n = 9;
x0 = [1;1;1;1;1;1;1;1;1];
N = 50;
%w为权重值
w = 1.25;
for i = 1:N
    [x, n] = sor(A,b,x0,w,i);
    fprintf('第%d次sor_%.2f迭代计算的结果:\n',n-1,w);
    disp(x);
end

% min = -20;   % 产生随机数的最小值
% max = 20;    % 产生随机数的最大值
% nrand = 20;
% Arand = round(min+(max-min)*rand(nrand,nrand));
% brand = round(min+(max-min)*rand(nrand,1));
% x0rand = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1];
% Nrand = 50;
% for i = 1:Nrand
%     [xrand, nrand] = sor(Arand,brand,x0rand,w,i);
%     fprintf('第%d次sor_%.2f迭代计算的结果:\n',nrand-1,w);
%     disp(xrand);
% end

function [x,n] = sor(A,b,x,w,it_max)
%  求线性方程组的sor(successive over-relaxation)迭代法，调用格式为
%  [x, n] = sor(A,b,x,w,it_max)
%  其中, A 为线性方程组的系数矩阵，w为权重值，b 为常数项
%  it_max 为最大迭代次数
%  x 为线性方程组的解，n-1为迭代次数

    D = diag(diag(A));%求A的对角矩阵
    L = tril(A,-1);%求A的下三角矩阵
    U = triu(A,1);%求A的上三角矩阵
    M = (D+w.*L)\((1-w).*D-w.*U);
    beta = w.*((D+w.*L)\b);

    n = 1;%迭代次数
    while n <= it_max
        x_after = diag(diag(M*x+beta));
        x = x_after;      %储存计算的值
        n = n+1;
    end
end
