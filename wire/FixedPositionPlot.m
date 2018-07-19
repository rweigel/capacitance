clear
load FixedPosition.mat

for k = 1:length(L)
    N(k) = length(L{k});
    L1(k) = L{k}(1);
    Lm(k) = L{k}(end/2);
end

figure(1);clf;grid on;
    semilogx(N,L1.*N,'b.','MarkerSize',20);
    hold on;grid on;
    semilogx(N,Lm.*N,'g.','MarkerSize',20);
    legend('\lambda_{end}/\lambda_{mid}','\lambda_{end}/\lambda_{mid} Bonifirm','Location','NorthWest');
    xlabel('N');

l0 = 0.5-0.152./log(N)-0.123./(log(N)).^2;
l1 = ( 0.0719+0.912./log(N)-0.874./(log(N)).^2 ).*log(N);
    
figure(2);clf;grid on;
    semilogx(N,L1./Lm,'k.','MarkerSize',30);
    hold on;
    semilogx(N,l1./l0,'r.','MarkerSize',30);
    grid on;
    hold on;
    legend('\lambda_{end}/\lambda_{mid}','\lambda_{end}/\lambda_{mid} Bonifirm','Location','NorthWest');
    xlabel('N');
    %legend('\lambda_{end}/\lambda_{mid}','Location','NorthWest');

figure(3);clf;grid on;
    loglog(N,L1./Lm,'k.','MarkerSize',30);
    hold on;
    loglog(N,l1./l0,'r.','MarkerSize',30);
    grid on;
    hold on;
    legend('\lambda_{end}/\lambda_{mid}','\lambda_{end}/\lambda{mid} Bonifirm','Location','NorthWest');
    %legend('\lambda_{end}/\lambda_{mid}','Location','NorthWest');


figure(4);clf;grid on;
    plot(N,L1./Lm,'k.','MarkerSize',30);
    hold on;
    plot(N,l1./l0,'r.','MarkerSize',30);
    grid on;
    hold on;
    legend('\lambda_{end}/\lambda_{mid}','\lambda_{end}/\lambda{mid} Bonifirm','Location','NorthWest');
    %legend('\lambda_{end}/\lambda_{mid}','Location','NorthWest');
    
break    
%saveplots('./figures/Capacitance_Force_Method_Density');

break

No = 7;
N = 10:100:10000;
x = log10(N);
x = x(1:end-1);
y = log10(diff(R1));
P = polyfit(x(No:end)',y(No:end)',1)
Ni = 1:1e5;
dsi = polyval(P,log10(Ni));

del = (-1-P(1));
Rinf = (10^P(2)/del)*(1/No^del)

plot(x,y,'.','MarkerSize',20);hold on;

plot(x,dsi,'-');hold on;
