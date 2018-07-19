function [Fss,x] = pushcode_err(p,N,lx)

xo = linspace(-1,1,N-1);
dx = (p(1)+log(1-xo.^2))/p(2);

dx(1) = interp1(xo(2:4),dx(2:4),xo(1),'pchip','extrap');
dx(end) = interp1(xo(end-3:end-1),dx(end-3:end-1),xo(end),'pchip','extrap');

x = [0,cumsum(dx)];
x = (x-min(x));
x = 2*x/max(x)-1;

if (nargin == 3)
    x = x+p(3)*lx;
end

%x = x+p(4)*sin(2*pi*[0:N-1]/N);

for i = 1:N
    F(i) = 0;
    for j = 1:N
        if i ~= j
            F(i) = F(i) + sign(x(i)-x(j))/(x(i)-x(j))^2;
        end
    end
end

% Boundary forces are not zero.
Fss = sqrt(sum(F(2:end-1).^2))/(length(F)-2);



