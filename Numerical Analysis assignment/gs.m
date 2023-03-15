%用G-S正交化的QR分解求解方程组Ax=b
function [x3,t] = gs(n)
tic;%计时开始
B=diag(ones(n,1))+diag(ones(n-1,1),-1);
A=B(1:n,1:n-1);%将A表示出来
C=eye(n);
b=zeros(n,1);
for k=1:n
    b=b+(((-1)^k)/n)*C(:,k);
end
b=b+2*ones(n,1)-C(:,1)-C(:,n);%将b表示出来
Q(:,1)=A(:,1)/norm(A(:,1),2);%第一个列向量单位化
for j=1:n-1 %一列一列地正交化
    y=A(:,j);
    for i=1:j-1
        R(i,j)=Q(:,i)'*A(:,j);
        y=y-R(i,j)*Q(:,i);
    end
    R(j,j)=norm(y,2);
    Q(:,j)=y/R(j,j);%将所求得正交基单位化
end
y=Q'*b;%求解Rx=y
for j=n-1:-1:2 %开始回代
    y(j)=y(j)/R(j,j);
    y(1:j-1)=y(1:j-1)-y(j)*R(1:j-1,j);
end
y(1)=y(1)/R(1,1);
x3=y;%得出解储存到x3
t=toc;%计时结束
end

