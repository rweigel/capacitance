clear;
load mat/pushcode_3495.mat

k = 1;
figure(1);clf;
for i = 100:length(X)
    if ~isempty(X{i})
        d     = diff(X{i});
        d1(k) = d(1);
        m     = ceil(length(d)/2);
        dm(k) = d(m);
        R(k) = dm(k)/d1(k);
        N(k) = i;
        Y{k} = X{i};
        %atan(sin(x)/(1+e^(-0.2*x)-cos(x)))
        %hold on;
        %I = find(X{i} < -0.99);
        I = find(X{i} < -1+dm(k)*2);
        f(k) = length(I)/N(k);
        %f(k) = d(m)/d(m-1);
        
        %loglog(N(k),sum(d(1:end/10)),'MarkerSize',30);drawnow;hold on;
        k = k+1;
    end
end

figure(1);clf;
loglog(N,d1,'.');
hold on;
loglog(N,dm,'.');
legend('\Delta_{end}','\Delta_{mid}');
grid on

figure(2);clf;
loglog(N(1:end-1),diff(log(d1)),'.');
hold on;
loglog(N(1:end-1),diff(log(dm)),'.');
legend('dlog(\Delta_{end})/dN','dlog(\Delta_{mid})/dN');
grid on

y1 = diff(log(d1));
x1 = diff(log(N));

ym = diff(log(dm));
xm = diff(log(N));

figure(3);clf;
loglog(N(1:end-2),-diff(y1)./diff(x1),'.');hold on;
loglog(N(1:end-2),-diff(ym)./diff(xm),'.');
legend('slope end','slope mid');
grid on

break

if (0)
I = [1:length(d1)];
c = 1e-3;
L(:,1) = 1./I;
L(:,2) = 1./sqrt(I);
L(:,3) = 1./log(I);
L(:,4) = exp(-c*I.^2);
L(:,5) = 1./I.^0.25;
L(:,6) = 1./(1+0.001*I.^2);

clf
plot(diff(d1),'LineWidth',2);hold on;
break


for i = 1:size(L,2),L(:,i) = L(:,i)/max(L(:,i));,end

clf
plot(d1/max(d1),'LineWidth',2);hold on;
plot(L,'LineWidth',2);
legend('\delta_1','1/i','1/sqrt(i)','1/log(i)','exp(-c*i^2)','1/i^{0.25}','1/(1+c*i^2)');

set(findall(gcf,'-property','FontSize'),'FontSize',16)

break


x = linspace(0,1,N(3000));
L(:,1) = (1/pi)*s./(s^2+x.^2);

s = 0.01;
L(:,2) = (1/s)*(1/sqrt(2*pi)).*exp(-x.^2/(2*s^2));

loglog(x,L/N(end));

hold on;
end

x = linspace(0,1,N(3000));
d = diff(X{3000});
I = find(x<0.25);
d = d(I);
ds = -d(I)-min(-d(I));
loglog(x(I),ds/ds(1),'.');

hold on;

x = linspace(0,1,N(300));
d = diff(X{300});
I = find(x<0.25);
d = d(I);
ds = -d(I)-min(-d(I));
loglog(x(I),0.7*ds/ds(1),'.')
