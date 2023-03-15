%用Householder正交化的QR分解求解方程组Ax=b
function [x4,t] = householder(n)
tic;%计时开始
B=diag(ones(n,1))+diag(ones(n-1,1),-1);
A=B(1:n,1:n-1);%将A表示出来
C=eye(n);
b=zeros(n,1);
for k=1:n
    b=b+(((-1)^k)/n)*C(:,k);
end
b=b+2*ones(n,1)-C(:,1)-C(:,n);%将b表示出来
z=b;
H=eye(n);
for i=1:n-1%计算householder变换矩阵
    [v,b]=house(A(i:n,i));
    G=eye(n+1-i)-b*(v')*v;
    P=eye(n);
    P(i:n,i:n)=G;
    H=P*H;
    A=P*A;
end
Q=H';
R=A;b=z;
y=Q'*b;%下面求解Rx=y
for j=n-1:-1:2
    y(j)=y(j)/R(j,j);
    y(1:j-1)=y(1:j-1)-y(j)*R(1:j-1,j);
end
y(1)=y(1)/R(1,1);
x4=y(1:n-1);%得出的解储存到x4
t=toc;%计时结束
end


