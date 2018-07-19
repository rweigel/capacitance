clear
load pushcode.mat
%'Display','iter',
options = optimset('TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',1000);
po = [9.75,1000,-1e-3];

n = 761;
dX = diff(X{n});
Ferr(X{n})

p = fminsearch(@(p) pushcode_err2(p,n),po,options);
[Fss,x3] = pushcode_err2(p,n);
Fss


dX = diff(X{n+1});
Ferr(X{n+1})

lx = X{n}-x3;
lxi = interp1(linspace(-1,1,n),lx,linspace(-1,1,n+1),'spline');

p = fminsearch(@(p) pushcode_err2(p,n+1,lxi),p,options);
[Fss,x3] = pushcode_err2(p,n+1,lxi);
Fss

f = @(a,x3,n) Ferr(x3+a*sin(2*pi*[0:n]/(n+1)));
a = fminsearch(@(a) f(a,x3,n),3e-5);
f(a,x3,n)
x4 = x3+a*sin(2*pi*[0:n]/(n+1));

figure(3);clf;hold on;
    %plot(dx,'b');
    plot(dX,'g','LineWidth',3);
    %plot(dx2,'k');
    plot(diff(x3),'m.','MarkerSize',10);

figure(4);clf;hold on;
   plot(x3-X{n+1},'k');
   plot(x4-X{n+1},'g')
   %plot(x5-X{n+1},'r')
   break 
figure(4);clf;hold on;
    plot(X{n},'g','LineWidth',3);
    plot(x,'m');
    
if (1)
    %options = optimset('TolFun',1e-9,'TolX',1e-9,'Display','iter','MaxFunEvals',50000);
    options = optimset('TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',50000);
    p = fminsearch(@(p) pushcode_err(p,n),po,options);
    [Fss,x] = pushcode_err(p,n);
    dx = diff(x);
    Ferr(x)
end

break

xo = linspace(-1,1,n);
dx = (p(1)+log(1-xo.^2)/p(2));
dx(1) = dx(2);
dx(end) = dx(end-1);
x = cumsum(dx);
x = (x-min(x));
x = 2*x/max(x)-1;

xo = linspace(-1,1,n-1);
D = diff(X{n});
px = polyfit(linspace(-1,1,length(D)),D,9);
tmpx = polyval(px,xo);
off = [0,cumsum(tmpx)];
xo = linspace(-1,1,n);
x2 = -1+2*off/off(end);
dx2 = diff(x2);
Ferr(x2)



