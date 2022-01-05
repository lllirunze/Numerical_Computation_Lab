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

% min = -20;   % 产生随机数的最小值
% max = 20;    % 产生随机数的最大值
% A = round(min+(max-min)*rand(nrand,nrand));
% b = round(min+(max-min)*rand(nrand,1));
  
N=length(b);    %解向量的维数
x=zeros(N,1);   %迭代近似向量
eps=0.0000001;  %精度
r=b-A*x;
d=r;
for k=0:N-1
    % fprintf('第%d次迭代：',k+1);
    a=(norm(r)^2)/(d'*A*d);
    x=x+a*d;
    rr=b-A*x;
    if (norm(rr)<=eps)||(k==N-1)
        break;
    end
    B=(norm(rr)^2)/(norm(r)^2);
    d=rr+B*d;
    r=rr;
end
disp(x)


