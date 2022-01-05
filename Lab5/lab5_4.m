clc
clear
close

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
x0 = [0;0;0;0;0;0;0;0;0];
[x,k] = max_speed_method(A,b,x0);
disp(x);

% min = -20;   % 产生随机数的最小值
% max = 20;    % 产生随机数的最大值
% nrand = 20;
% Arand = round(min+(max-min)*rand(nrand,nrand));
% brand = round(min+(max-min)*rand(nrand,1));
% x0rand = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
% [xrand,krand] = max_speed_method(Arand,brand,x0rand);

function [x,k] = max_speed_method(A,b,x0)%最速下降法
    k=1;%迭代次数
    N=1000;%限制最大迭代次数
    x=x0;
    while k<N%迭代过程
        r=b-A*x;%步长
        l=dot(r,r)/dot(A*r,r);%方向
        x=x+l*r;
        if norm(x-x0)<10^(-4)%精度
            break;
        else
            k=k+1;
        end
        x0=x;
    end
end
