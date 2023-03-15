%用法方程组的共轭梯度法求解方程组Ax=b
function [x2,t] = fagonge(n)
tic;%计时开始
B=diag(ones(n,1))+diag(ones(n-1,1),-1);
A=B(1:n,1:n-1);%将A表示出来
C=eye(n);
b=zeros(n,1);
for k=1:n
    b=b+(((-1)^k)/n)*C(:,k);
end
b=b+2*ones(n,1)-C(:,1)-C(:,n);%将b表示出来
b=A'*b;A=A'*A;
n=n-1;%将n换一种表示方法，使后面的语句更简洁一点
x2=zeros(n,1);%选取初值为0向量
r=b-A*x2;k=0;m=r;%计算初始残差，准备开始迭代
while norm(r,2)>0.00001%确定残差范数精度
    k=k+1;
    if k==1%第一步最速下降
        p=r;
    else
        l=r;r=m;m=l;%然后共轭梯度
        w=m'*m/(r'*r);
        p=m+w*p;
    end
    v=m'*m/(p'*A*p);
    x2=x2+v*p;%迭代若干次后得出所求解
    r=m-v*A*p;
end
t=toc;%计时结束
end

