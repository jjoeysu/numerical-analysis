%Householder变换过程
function [v,b]=house(x)
n =length(x);
i=norm(x,inf);
x =x/i;
d =x(2:n)'*x(2:n);
v(2:n) =x(2:n);
if d==0 
    b =0;
else
    a =sqrt(x(1)^2+d);
    if x(1)<=0
        v(1) =x(1)-a;
    else
        v(1)=-d/(x(1)+a);
    end
    b=2*v(1)^2/(d+v(1)^2);
    v=v/v(1);
end
end
