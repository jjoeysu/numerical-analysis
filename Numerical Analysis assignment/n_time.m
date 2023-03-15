%下面是各个方法所用计算时间随n的变化曲线
N=100;%确定n的最大值N

%法方程组的cholesky分解
%经计算，其计算量约为（10/3）*n^3
for i=2:N
n1(i)=i;
[x1,t] = facholesky(i);
time1(i)=t;
end
subplot(2,2,1);
plot(n1,time1);
xlabel('n'),ylabel('所用时间 t / s');
title('法方程组的cholesky分解');

%法方程组的共轭梯度法(初值为0向量)
for i=2:N
n2(i)=i;
[x2,t] = fagonge(i);
time2(i)=t;
end
subplot(2,2,2);
plot(n2,time2);
xlabel('n'),ylabel('所用时间 t / s');
title('法方程组的共轭梯度法');

%G-S正交化的QR分解
%经计算，其计算量约为6*n^3
for i=2:N
n3(i)=i;
[x3,t] = gs(i);
time3(i)=t;
end
subplot(2,2,3);
plot(n3,time3);
xlabel('n'),ylabel('所用时间 t / s');
title('G-S正交化的QR分解');

%Householder变换的QR分解
%经计算，其计算量约为（4/3）*n^3
for i=2:N
n4(i)=i;
[x4,t] = householder(i);
time4(i)=t;
end
subplot(2,2,4);
plot(n4,time4);
xlabel('n'),ylabel('所用时间 t / s');
title('Householder变换的QR分解');