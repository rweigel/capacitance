clear;
load mat/pushcode_3495.mat

%load mat/pushcodec.mat
%X = DataC.X;

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
semilogx(N,f);grid on;
title(f(end));
break
No = 7;
x = log10(N(1:end-1));
y = log10(diff(R));
P = polyfit(x(No:end)',y(No:end)',1)
Ni = N(No:end):1e4;
dsi = polyval(P,log10(Ni));

del = (-1-P(1));
Rinf = (10^P(2)/del)*(1/No^del);
break
if (0)
figure(2);clf;hold on;
    plot(log10(N),log(d1),'r.','LineWidth',3);grid on;
    plot(log10(N),log(dm),'g.','LineWidth',3);grid on;
    plot(log10(N),log(1./R),'b.','LineWidth',3);grid on;    
    xlabel('log_{10}(N)');
    legend('\Delta_{end}','\Delta_{mid}','R\equiv\Delta_{mid}/\Delta_{end}','Location','SouthWest');
end

if (1)
figure(2);clf;hold on;
    plot(log10(N),1./(N.*d1),'r.','MarkerSize',20);grid on;
    plot(log10(N),1./(N.*dm),'g.','MarkerSize',20);
    %plot(log10(N),dm./d1,'b','LineWidth',3);
    xlabel('log_{10}(N)');
    legend('\lambda_{end}','\lambda_{mid}','Location','NorthWest');
    %legend('\lambda_{end}','\lambda_{mid}','R\equiv\lambda_{end}/\lambda_{mid}','Location','NorthWest');
end

figure(3);clf;hold on;
    plot(log10(N(1:end-1)),log10(diff(R)),'k.','MarkerSize',10);
    plot(log10(Ni),dsi,'g','LineWidth',1);
    xlabel('log_{10}(N)');
    legend('dR/dN','dR/dN interpolated');
    title(sprintf('slope = %.3f; R_{\\infty} = %.2f',P(1),Rinf));
    grid on;

break    
figure(4);clf;
k = 1;
for i = [5:length(Y)]
    d = diff(Y{i});
    D = d(1:end-1)./d(2:end);
    N = 1:length(D);
    loglog((N(1:end/2))/N(end/2),log(1./D(1:end/2)),'.');
    hold on;
    %L{k} = sprintf('N = %d',length(X{end})-i);
    Nm(k) = N(end/2);
    D1(k) = D(1);
    Dm(k) = 1./D(end/2);
    Dmm(k) = 1./D(round(end/2-1));
    polyfit(log(N(end/2-1:end/2)),log(log(1./D(end/2-1:end/2))),1)
    k = k+1;
end
xlabel('i');
title('log_{10}(\Delta_{i+1}/\Delta_{i})')
%legend(L{:})
grid on;
polyfit(log(Nm),log(log(Dm)),1)
polyfit(log(Nm),log(log(Dmm)),1)

% Seems to follow 1/(1+x^a)
break
n = N(1:end/2);

figure(7);clf;
loglog(log(1./D(1:end/2)),N(1:end/2),'k.');
hold on;
loglog(1./(10+n),n)
break
    %plot(log(log10(N)),1./ds);


figure(3);clf;hold on;
    plot((log10(N)),ds);
    plot(log10(N),(0.9*l1./l0));
    hold on;grid on;
    %plot(log(log10(N)),1./ds);

%plot(5*N/(max(N)),exp(-(N/max(N)).^0.33))
%plot(1.5+N/(max(N)),exp(-(N/max(N)).^0.33))
xlabel('log(N)')
ylabel('\lambda_{0}/\lambda_{1}');
exp(-1)

%for i = 3:length(X)
%    FSS(i) = Ferr(X{i});
%end
if (0)
    xg = linspace(-1,1,1e7);
    xg = linspace(X{500}(200),X{500}(201),1e4);
    V = potential(xg,X{500});
    clf
    semilogy(xg,V);hold on
    %plot(X{500},3600*ones(500,1),'.')
    set(gca,'YLim',[3e3,4e3]);
    clf
    I = find(xg>0);
    loglog(V(1:end/2),'.')
    grid on
    plot(Vo+Vo2-V)
end


k = 1;
for i = 12:length(X)
    if ~isempty(X{i})
        D = diff(X{i});
        tmp = polyfit(linspace(-1,1,length(D)),D,13);
        tmp = tmp(2:2:end);
        P(k,1:length(tmp)) = tmp;
        k = k+1;
    end
end
figure(1);
    plot((P),'LineWidth',2);
    grid on;
    for i = 1:size(P,2)
        text(k+1,P(end,i),sprintf('c_{%2d}',2*i));
    end
    


%p = [0.91,5,1,1,0.1];
p = [0.91,5,1];
%p = [0.91,5,1];
break
for i = 210:210%length(X)
    if ~isempty(X{i})
        d = diff(X{i});
        %options = optimset('TolFun',1e-9,'TolX',1e-9,'Display','iter','MaxFunEvals',50000);
        options = optimset('TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',50000);
        p = fminsearch(@(p) xerrfn(p,i),p,options);
        [Fmse,x] = xerrfn(p,i);
        fprintf('%03d %g\n',i,Fmse);
        fmse(i) = Fmse;
        figure(3);
            plot(diff(X{i})-diff(x))
            plot(X{i}(2:end-1)-x(2:end-1));
            title(sprintf('n=%d',i));
        y = (x(2:end-1)-X{i}(2:end-1));
        y(1)=0;
        y(end)=0;
        I = fft(y).*conj(fft(y));
        figure(4);
            loglog(I(2:end/2-1),'*')
            title(sprintf('n=%d',i));
        
        %plot(linspace(-1,1,length(x)-1),diff(x),'k');
        cc=corrcoef(d,diff(x));
        CC(i)=cc(2);
        drawnow
        k = k+1;
    end
end

    
break

figure(3);clf;hold on;grid on;
xlabel('x');
ylabel('p-x');
title('Deviation from linear')
for n = 3:length(X)
    plot(X{n},Y{n},'-');
end

figure(4);clf;hold on;grid on;box on;
xlabel('charge #/N')
ylabel('x')
for n = 3:length(X)
    I = [1:n];
    xgrid = linspace(-1,1,n);
    plot(I,Y{n}+xgrid,'-');
end
break

figure(5);clf;hold on;grid on;
xlabel('charge #/N');
ylabel('(x_{i+1}-x_{i})/dx');
%axis([0,1,0.5,1]);
%title('')
for n = 3:length(X)
    I = [0:n-1];
    xgrid = linspace(-1,1,n);
    D = diff(Y{n}+xgrid);
    plot((I(2:end)-1/2)/n,D/n)
    %pause(0.1);
end

figure(6);clf;hold on;grid on;
for n = 30:length(X)
    I = [0:n-1];
    xgrid = linspace(-1,1,n);
    D = diff(Y{n}+xgrid);
    x = linspace(-1,1,length(D));
    dx = x(2)-x(1);
    x = linspace(-1+dx/2,1-dx/2,length(D));
    xlabel('x');
    %ylabel('$\frac{\Delta}{\Delta_1}$','Interpreter','Latex','Rotation',0);
    ylabel('${\Delta}/{\Delta_1}$','Interpreter','Latex');
    title('\Delta = x_{i+1}-x_{i}','FontWeight','normal');
    P = polyfit(x,D,20);
    Df = polyval(P,x);
    %mD = max(D);
    mD = D(1);
    %plot(x,D/mD,'k','LineWidth',2);
    %plot(x,Df/mD,'w');
    plot(abs(P(P>0)))
    pause(0.5);
end

