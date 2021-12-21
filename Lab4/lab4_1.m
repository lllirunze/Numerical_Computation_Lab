clc
clear
close

syms x
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
% 列主元素消去法
[x] = Column_Pivot_Elimination(A,b);

syms xrand
n = 20;
min = -20;   % 产生随机数的最小值
max = 20;    % 产生随机数的最大值
Arand = round(min+(max-min)*rand(n,n));
brand = round(min+(max-min)*rand(n,1));
[xrand] = Column_Pivot_Elimination(Arand,brand);


function [x] = Column_Pivot_Elimination(A,b)
    n=length(A);
    [x]=zeros(n,1);
    
    for k=1:n
        maxindex=k;
        for u=k+1:n
            if(abs(A(u,k))>abs(A(k,k)))
                maxindex=u;
            end
        end
        
        temp=A(maxindex,:);
        A(maxindex,:)=A(k,:);
        A(k,:)=temp;
        
        bt=b(maxindex);
        b(maxindex)=b(k);
        b(k)=bt;
        
        r=A(1:k,1:k);
        Det=det(r);
        if Det==0
            error ('Algorithm can not solve this matrix.');
        end
        
        for i=k+1:n
            mi=A(i,k)/A(k,k);
            b(i)=b(i)-mi*b(k);
            for j=k+1:n
                A(i,j)=A(i,j)-mi*A(k,j);
            end
        end
    end
    x(n)=b(n)/A(n,n);
    for i=n-1:-1:1
        x(i)=(b(i)-sum(A(i,i+1:n)*x(i+1:n)))/A(i,i);
    end
end

