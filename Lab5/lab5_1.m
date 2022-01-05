clc
clear
close

% 输入量：
% A: 线性方程组的系数矩阵（n*n，非奇异）
% b: 方程组右边的常数项列向量
% n: 方程组维数
% x0: 初始值
% tol: 精度上限值
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
x0 = [0;0;0;0;0;0;0;0;0];
tol = 1e-6;
N = 50;
x = jacobi_fun(A,b,n,x0,tol,N);

min = -20;   % 产生随机数的最小值
max = 20;    % 产生随机数的最大值
nrand = 20;
Arand = round(min+(max-min)*rand(nrand,nrand));
brand = round(min+(max-min)*rand(nrand,1));
x0rand = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
tolrand = 1e-6;
Nrand = 50;
xrand = jacobi_fun(Arand,brand,nrand,x0rand,tolrand,Nrand);


function x = jacobi_fun(A,b,n,x0,tol,N)
x = zeros(n,1);      % 给x赋值
k = 0;
while k<N
    for i=1:n
     x(i)=(b(i)-A(i,[1:i-1,i+1:n])*x0([1:i-1,i+1:n]))/A(i,i);
    end
    if norm(x-x0)<tol
        break;
    end
    x0=x;
    k=k+1;

    disp(['when k=',num2str(k)])
    disp('x=');
    disp(x);                       %输出计算的中间结果
end

if k==N
    disp('迭代次数已到达上限!');
end
disp(['迭代次数 k=',num2str(k)])

end
