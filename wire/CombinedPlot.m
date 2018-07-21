clear;

%Data = load('mat/FixedPosition_100.mat');
%Lambda = Data.Lambda;

Data = load('mat/FixedPosition_32768.mat');
Lambda = Data.L;

for k = 1:length(Lambda)
    Nfp(k) = length(Lambda{k});
    % Lambda{k} is really charge on each patch. Total charge is always 1,
    % length of interval is 2.
    % Lambda = q_i/dx = q_i/ (2/Nfp);
    L1fp(k) = Lambda{k}(1)*Nfp(k)/2;
    D1fp(k) = (1./Nfp(k)).*(2/Nfp(k))./Lambda{k}(1); % dx = dq/lambda
    Lmfp(k) = Lambda{k}(end/2)*Nfp(k)/2;
    Dmfp(k) = (1./Nfp(k)).*(2/Nfp(k))./Lambda{k}(end/2); % dx = dq/lambda
end

%load mat/pushcodec.mat
%X = DataC.X;

Data = load('mat/pushcode_3495.mat');
X = Data.X;

k = 1;
for i = 10:length(X)
    if ~isempty(X{i})
        Xfc{k} = X{i};
        d      = diff(X{i});
        Nfc(k) = i;
        % Assume total charge is 1.
        % Lamda = dq/dx = (1/N)/(x_{i+1}-x_{i})
        L1fc(k) = 1./(Nfc(k)*d(1));
        D1fc(k) = d(1);
        m       = ceil(length(d)/2);
        Lmfc(k) = 1./(Nfc(k)*d(m));
        Dmfc(k) = d(m);
        k = k+1;
    end
end

figure(1);clf;grid on;hold on;
    for i = [1000:100:length(Xfc)]
        d = diff(Xfc{i});
        x = d(1)/2+Xfc{i}(1:end-1); % Center point between each charge.
        plot(x,d);
    end
    legend off
    title(sprintf('N = [1000:100:%d]',length(Xfc)));
    xlabel('x');
    set(gca,'xlim',[-1.02,1.02])
    ylabel('\Delta (Fixed Charge)');
    set(findall(gcf,'-property','FontSize'),'FontSize',16)

figure(2);clf;grid on;hold on;
    for i = [1000:100:length(Xfc)]
        d = diff(Xfc{i});
        x = d(1)/2+Xfc{i}(1:end-1); % Center point between each charge.
        plot(x,d,'Marker','.','LineStyle','-');
    end
    legend off
    title(sprintf('N = [1000:100:%d]',length(Xfc)));
    xlabel('x');
    set(gca,'xlim',[-1.001,-0.99])
    ylabel('\Delta (Fixed Charge)');
    set(findall(gcf,'-property','FontSize'),'FontSize',16)

figure(3);clf;grid on;
    semilogx(Nfc,D1fc,'k.','MarkerSize',20);
    hold on;grid on;
    semilogx(Nfp,D1fp/0.873,'r.','MarkerSize',20);
    semilogx(Nfc,Dmfc,'k.','MarkerSize',10);
    semilogx(Nfp,Dmfp,'r.','MarkerSize',10);
    legend('\Delta_{end} (Fixed Charge)','\Delta_{end} (Fixed Position)','\Delta_{mid} (Fixed Charge)','\Delta_{mid} (Fixed Position)','Location','NorthEast');
    xlabel('N');
    %axis tight;
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    
figure(4);clf;grid on;
    loglog(Nfc,D1fc,'k.','MarkerSize',20);
    hold on;grid on;
    loglog(Nfp,D1fp/0.873,'r.','MarkerSize',20);
    loglog(Nfc,Dmfc,'k.','MarkerSize',10);
    loglog(Nfp,Dmfp,'r.','MarkerSize',10);
    legend('\Delta_{end} (Fixed Charge)','\Delta_{end} (Fixed Position)','\Delta_{mid} (Fixed Charge)','\Delta_{mid} (Fixed Position)','Location','NorthEast');
    xlabel('N');
    axis tight;
    set(findall(gcf,'-property','FontSize'),'FontSize',16)

figure(5);clf;grid on;
    semilogx(Nfc,D1fc./Dmfc,'k.','MarkerSize',20);
    hold on;grid on;
    semilogx(Nfp,(1/0.873)*D1fp./Dmfp,'r.','MarkerSize',20);
    legend('\Delta_{end}/\Delta_{mid} (Fixed Charge)','(1/0.873)\Delta_{end}/\Delta_{mid} (Fixed Position)','Location','NorthEast');
    xlabel('N');
    %axis tight;
    set(gca,'ylim',[0.3,0.75]);
    set(findall(gcf,'-property','FontSize'),'FontSize',16)

figure(6);clf;grid on;
    loglog(Nfc,D1fc./Dmfc,'k.','MarkerSize',20);
    hold on;grid on;
    loglog(Nfp,(1/0.873)*D1fp./Dmfp,'r.','MarkerSize',20);
    legend('\Delta_{end}/\Delta_{mid} (Fixed Charge)','(1/0.873)\Delta_{end}/\Delta_{mid} (Fixed Position)','Location','NorthEast');
    xlabel('N');
    %axis tight;
    set(gca,'ylim',[0.3,0.75]);
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    
figure(7);clf;grid on;
    semilogx(Nfc,L1fc,'k.','MarkerSize',20);
    hold on;grid on;
    semilogx(Nfp,L1fp*0.873,'r.','MarkerSize',20);
    semilogx(Nfc,Lmfc,'k.','MarkerSize',10);
    semilogx(Nfp,Lmfp,'r.','MarkerSize',10);
    legend('\lambda_{end} (Fixed Charge)','0.873*\lambda_{end} (Fixed Position)','\lambda_{mid} (Fixed Charge)','\lambda_{mid} (Fixed Position)','Location','NorthWest');
    xlabel('N');
    axis tight;
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    
figure(8);clf;grid on;
    loglog(Nfc,L1fc,'k.','MarkerSize',20);
    hold on;grid on;
    loglog(Nfp,L1fp*0.873,'r.','MarkerSize',20);
    loglog(Nfc,Lmfc,'k.','MarkerSize',10);
    loglog(Nfp,Lmfp,'r.','MarkerSize',10);
    legend('\lambda_{end} (Fixed Charge)','0.873*\lambda_{end} (Fixed Position)','\lambda_{mid} (Fixed Charge)','\lambda_{mid} (Fixed Position)','Location','NorthWest');
    xlabel('N');
    axis tight;
    set(findall(gcf,'-property','FontSize'),'FontSize',16)

figure(9);clf;grid on;
    semilogx(Nfc,L1fc./Lmfc,'k.','MarkerSize',20);
    hold on;grid on;
    semilogx(Nfp,(0.873)*L1fp./Lmfp,'r.','MarkerSize',20);
    legend('\lambda_{end}/\lambda_{mid} (Fixed Charge)','(1/0.873)\lambda_{end}/\lambda_{mid} (Fixed Position)','Location','NorthWest');
    xlabel('N');
    %axis tight;
    set(gca,'ylim',[1.5,3.0]);
    set(findall(gcf,'-property','FontSize'),'FontSize',16)

figure(10);clf;grid on;
    loglog(Nfc,L1fc./Lmfc,'k.','MarkerSize',20);
    hold on;grid on;
    loglog(Nfp,(0.873)*L1fp./Lmfp,'r.','MarkerSize',20);
    legend('\lambda_{end}/\lambda_{mid} (Fixed Charge)','(1/0.873)\lambda_{end}/\lambda_{mid} (Fixed Position)','Location','NorthWest');
    xlabel('N');
    %axis tight;
    set(gca,'ylim',[1.5,3.0]);
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
 
break

% From Bonfim and Griffiths 2000
l0 = 0.5-0.152./log(Nfp)-0.123./(log(Nfp)).^2;
l1 = ( 0.0719+0.912./log(Nfp)-0.874./(log(Nfp)).^2 ).*log(Nfp);

options = optimset('TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',50000);
po = [0.07,1.0,-0.874];
po = [0.07,1.0];
p = fminsearch(@(p) error_fit(p,Nfp,L1fp),po,options);


L1fp_fit = p(1) + p(2).*log(Nfp);

figure(3);clf;grid on;
    semilogx(Nfp,0.873*L1fp,'b.','MarkerSize',30);
    hold on;grid on;
    semilogx(Nfc,L1fc,'r.','MarkerSize',20);
    semilogx(Nfp,0.84*l1,'g.','MarkerSize',20);
    semilogx(Nfp,0.873*L1fp_fit,'k-');
    legend('0.873\lambda_{end} (Fixed Position)','\lambda_{end} (Fixed Charge)','0.8\lambda_{end} (Bonifirm Fit)','0.873\lambda_{end} (Fixed Position Fit for N > 100)','Location','NorthWest');
    xlabel('N');
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
    
dL1fp = diff(L1fp)./diff(Nfp); 
dL1fc = diff(L1fc)./diff(Nfc);

no = 1000;
nf = 3000;
x = log10(Nfc(1:end-1));
y = log10(dL1fc);
P = polyfit(x(no:nf)',y(no:nf)',1)
dL1fc_fit = polyval(P,log10(Nfc));

% Model is dL1/dN = const*N^alpha = P(1)*N^P(2)
% Integral is L1 = [const/(alpha+1)]*N^(alpha+1)|_no^nf
% L1(inf) = -const/(alpha+1)*N^(alpha+1);

L1fc_inf = (-P(1)/(P(2)+1))*no^(P(2)+1)

figure(3);clf;grid on;
    loglog(Nfp(2:end),dL1fp,'b.','MarkerSize',30);
    hold on;grid on;
    loglog(Nfc(2:end),dL1fc,'r.','MarkerSize',20);
    loglog(Nfc,10.^dL1fc_fit,'g-','MarkerSize',30);
    legend('d\lambda_{end}/dN (Fixed Position)','d\lambda_{end}/dN (Fixed Charge)','Location','NorthWest');
    xlabel('N');
    set(findall(gcf,'-property','FontSize'),'FontSize',16)
   

    