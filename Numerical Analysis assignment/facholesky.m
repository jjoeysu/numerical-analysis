%用法方程组的cholesky分解求解方程组Ax=b
function [x1,t] = facholesky(n)
tic;%计时开始
B=diag(ones(n,1))+diag(ones(n-1,1),-1);
A=B(1:n,1:n-1);%将A表示出来
C=eye(n);
b=zeros(n,1);
for k=1:n
    b=b+(((-1)^k)/n)*C(:,k);
end
b=b+2*ones(n,1)-C(:,1)-C(:,n);%将b表示出来
b=A'*b;A=A'*A;%转换A和b,便于分解
%开始平方根法计算cholesky分解
for i=1:n-1
    A(i,i)=sqrt(A(i,i));
    A(i+1:n-1,i)=A(i+1:n-1,i)/A(i,i);
    for j=i+1:n-1
        A(j:n-1,j)=A(j:n-1,j)-A(j:n-1,i)*A(j,i);
    end
end
%开始前代
for l=1:n-2
    b(l)=b(l)/A(l,l);
    b(l+1:n-1)=b(l+1:n-1)-b(l)*A(l+1:n-1,l);
end
b(n-1)=b(n-1)/A(n-1,n-1);
U=A';
%开始回代
for m=n-1:-1:2
    b(m)=b(m)/U(m,m);
    b(1:m-1)=b(1:m-1)-b(m)*U(1:m-1,m);
end
b(1)=b(1)/U(1,1);
x1=b;
t=toc;%计时结束
end

