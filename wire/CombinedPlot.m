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
    Lfp{k}   = Lambda{k}*Nfp(k)/2;
    Nfp_a(k) = length(find(Lfp{k}>0.5));
    Nfp_b(k) = length(find(Lfp{k}<=0.5));
    Ffp_a(k) = Nfp_a(k)/Nfp(k); % Fraction above
    Ffp_b(k) = Nfp_b(k)/Nfp(k); % Fraction below
    Afp_a(k) = (1/2)*sum(Lfp{k}(Lfp{k}>0.5))/Nfp_a(k); % Area above
    Afp_b(k) = (1/2)*sum(Lfp{k}(Lfp{k}<=0.5))/Nfp_b(k); % Area above
    L1fp(k) = Lambda{k}(1)*Nfp(k)/2;
    D1fp(k) = (1./Nfp(k)).*(2/Nfp(k))./Lambda{k}(1); % dx = dq/lambda
    L2fp(k) = Lambda{k}(2)*Nfp(k)/2;
    D2fp(k) = (1./Nfp(k)).*(2/Nfp(k))./Lambda{k}(2); % dx = dq/lambda
    Lmfp(k) = Lambda{k}(end/2)*Nfp(k)/2;
    Dmfp(k) = (1./Nfp(k)).*(2/Nfp(k))./Lambda{k}(end/2); % dx = dq/lambda
    Lm2fp(k) = Lambda{k}(end/2-1)*Nfp(k)/2;
    Dm2fp(k) = (1./Nfp(k)).*(2/Nfp(k))./Lambda{k}(end/2-1); % dx = dq/lambda

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
        Lfc{k}   = 1./(Nfc(k)*d);
        xfc(k)   = X{i}(find(Lfc{k}<=0.5,1)); % First location where Lambda <= 0.5.
        Nfc_a(k) = length(find(Lfc{k}>0.5));
        Nfc_b(k) = length(find(Lfc{k}<=0.5));
        Ffc_a(k) = Nfc_a(k)/Nfc(k); % Fraction above
        Ffc_b(k) = Nfc_b(k)/Nfc(k); % Fraction below
        Afc_a(k) = (1/2)*sum(Lfc{k}(Lfc{k}>0.5))/Nfc_a(k); % Area above
        Afc_b(k) = (1/2)*sum(Lfc{k}(Lfc{k}<=0.5))/Nfc_b(k); % Area above
        L1fc(k) = 1./(Nfc(k)*d(1));
        D1fc(k) = d(1);
        L2fc(k) = 1./(Nfc(k)*d(2));
        D2fc(k) = d(2);
        m       = ceil(length(d)/2);
        Lmfc(k) = 1./(Nfc(k)*d(m));
        Dmfc(k) = d(m);
        Lm2fc(k) = 1./(Nfc(k)*d(m-1));
        Dm2fc(k) = d(m-1);
        k = k+1;
    end
end

fn = 0;

I = [100:400:length(Xfc)];
Is = sprintf('%d,',I);
Is = sprintf('N = %s',Is(1:end-1));

fn=fn+1;figure(fn);clf;grid on;hold on;
    for i = I
        d = Lfc{i};
        dx = diff(Xfc{i});
        x = dx/2+Xfc{i}(1:end-1); % Center point between each charge.
        plot(x,d);
    end
    legend off
    title(Is);
    xlabel('x');
    set(gca,'xlim',[-1.02,1.02])
    ylabel('$\lambda$ (Fixed Charge)');
    setfonts

fn=fn+1;figure(fn);clf;grid on;hold on;
    for i = I
        d = Lfc{i};
        dx = diff(Xfc{i});
        x = dx/2+Xfc{i}(1:end-1); % Center point between each charge.
        plot(x,d,'Marker','.','LineStyle','-');
    end
    legend off
    title(Is);
    xlabel('x');
    set(gca,'xlim',[-1.001,-0.99])
    ylabel('$\lambda$ (Fixed Charge)');
    setfonts
 
fn=fn+1;figure(fn);clf;grid on;hold on;
    for i = I
        d = diff(Xfc{i});
        x = d(1)/2+Xfc{i}(1:end-1); % Center point between each charge.
        plot(x,d);
    end
    legend off
    title(Is);
    xlabel('x');
    set(gca,'xlim',[-1.02,1.02])
    ylabel('$\Delta$ (Fixed Charge)');
    setfonts

fn=fn+1;figure(fn);clf;grid on;hold on;
    for i = I
        d = diff(Xfc{i});
        x = d(1)/2+Xfc{i}(1:end-1); % Center point between each charge.
        plot(x,d,'Marker','.','LineStyle','-');
    end
    legend off
    title(Is);
    xlabel('x');
    set(gca,'xlim',[-1.001,-0.99])
    ylabel('$\Delta$ (Fixed Charge)');
    setfonts

fn=fn+1;figure(fn);clf;
    semilogx(Nfc,D1fc,'k.','MarkerSize',20);
    hold on;grid on;
    semilogx(Nfp,D1fp/0.873,'r.','MarkerSize',20);
    semilogx(Nfc,Dmfc,'k.','MarkerSize',10);
    semilogx(Nfp,Dmfp,'r.','MarkerSize',10);
    legend('$\Delta_{0}$ (Fixed Charge)','$\Delta_{0}$ (Fixed Position)','$\Delta_{N/2}$ (Fixed Charge)','$\Delta_{N/2}$ (Fixed Position)','Location','NorthEast');
    xlabel('N');
    %axis tight;
    setfonts
    
fn=fn+1;figure(fn);clf;
    loglog(Nfc,D1fc,'k.','MarkerSize',20);
    hold on;grid on;
    loglog(Nfp,D1fp/0.873,'r.','MarkerSize',20);
    loglog(Nfc,Dmfc,'k.','MarkerSize',10);
    loglog(Nfp,Dmfp,'r.','MarkerSize',10);
    legend('$\Delta_{0}$ (Fixed Charge)','$\Delta_{0}$ (Fixed Position)','$\Delta_{N/2}$ (Fixed Charge)','$\Delta_{N/2}$ (Fixed Position)','Location','NorthEast');
    xlabel('N');
    axis tight;
    setfonts

fn=fn+1;figure(fn);clf;
    semilogx(Nfc,D1fc./Dmfc,'k.','MarkerSize',20);
    hold on;grid on;
    semilogx(Nfp,(1/0.873)*D1fp./Dmfp,'r.','MarkerSize',20);
    legend('$\Delta_{0}/\Delta_{N/2}$ (Fixed Charge)','$(1/0.873)\Delta_{0}/\Delta_{N/2}$ (Fixed Position)','Location','NorthEast');
    xlabel('N');
    set(gca,'ylim',[0.3,0.80]);
    setfonts

fn=fn+1;figure(fn);clf;
    semilogx(Nfc,D1fc./D2fc,'r.','MarkerSize',20);
    hold on;grid on;
    semilogx(Nfp,1.1501*D1fp./D2fp,'b.','MarkerSize',20);
    legend('$\Delta_{0}/\Delta_{1}$ (Fixed Charge)','$(1.1501)\Delta_{0}/\Delta_{1}$ (Fixed Position)','Location','NorthEast');
    xlabel('N');
    %set(gca,'ylim',[0.3,0.80]);
    setfonts
    
fn=fn+1;figure(fn);clf;
    loglog(Nfc,D1fc./Dmfc,'k.','MarkerSize',20);
    hold on;grid on;
    loglog(Nfp,(1/0.873)*D1fp./Dmfp,'r.','MarkerSize',20);
    legend('$\Delta_{0}/\Delta_{N/2}$ (Fixed Charge)','$(1/0.873)\Delta_{0}/\Delta_{N/2}$ (Fixed Position)','Location','NorthEast');
    xlabel('N');
    set(gca,'ylim',[0.3,0.80]);
    setfonts
    
fn=fn+1;figure(fn);clf;
    loglog(Nfc,D2fc./D1fc,'k.','MarkerSize',20);
    hold on;grid on;
    loglog(Nfp,(1/1.15)*D2fp./D1fp,'r.','MarkerSize',20);
    legend('$\Delta_{0}/\Delta_{1}$ (Fixed Charge)','$1.15\Delta_{0}/\Delta_{1}$ (Fixed Position)','Location','NorthEast');
    xlabel('N');
    %set(gca,'ylim',[0.995,1.01]);
    setfonts
    
fn=fn+1;figure(fn);clf;
    loglog(Nfc,1-Dm2fc./Dmfc,'k.','MarkerSize',20);
    hold on;grid on;
    loglog(Nfp,1-Dm2fp./Dmfp,'r.','MarkerSize',20);
    legend('$1-\Delta_{N/2-1}/\Delta_{N/2}$ (Fixed Charge)','$1-\Delta_{N/2-1}/\Delta_{N/2}$ (Fixed Position)','Location','NorthEast');
    xlabel('N');
    %set(gca,'ylim',[0.995,1.01]);
    setfonts
    
fn=fn+1;figure(fn);clf;
    semilogx(Nfc,L1fc,'k.','MarkerSize',20);
    hold on;grid on;
    semilogx(Nfp,L1fp*0.873,'r.','MarkerSize',20);
    semilogx(Nfc,Lmfc,'k.','MarkerSize',10);
    semilogx(Nfp,Lmfp,'r.','MarkerSize',10);
    legend('$\lambda_{0}$ (Fixed Charge)','$0.873\lambda_{0}$ (Fixed Position)','$\lambda_{N/2}$ (Fixed Charge)','$\lambda_{N/2}$ (Fixed Position)','Location','NorthWest');
    xlabel('N');
    set(gca,'ylim',[0.38,1.4]);
    setfonts
    
fn=fn+1;figure(fn);clf;
    loglog(Nfc,L1fc,'k.','MarkerSize',20);
    hold on;grid on;
    loglog(Nfp,L1fp*0.873,'r.','MarkerSize',20);
    loglog(Nfc,Lmfc,'k.','MarkerSize',10);
    loglog(Nfp,Lmfp,'r.','MarkerSize',10);
    legend('$\lambda_{0}$ (Fixed Charge)','$0.873\lambda_{0}$ (Fixed Position)','$\lambda_{N/2}$ (Fixed Charge)','$\lambda_{N/2}$ (Fixed Position)','Location','NorthWest');
    xlabel('N');
    set(gca,'ylim',[0.38,1.4]);
    setfonts

fn=fn+1;figure(fn);clf;
    semilogx(Nfc,L1fc./Lmfc,'k.','MarkerSize',20);
    hold on;grid on;
    semilogx(Nfp,(0.873)*L1fp./Lmfp,'r.','MarkerSize',20);
    legend('$\lambda_{0}/\lambda_{N/2}$ (Fixed Charge)','$(1/0.873)\lambda_{0}/\lambda_{N/2}$ (Fixed Position)','Location','NorthWest');
    xlabel('N');
    %axis tight;
    set(gca,'ylim',[1.2,3.0]);
    setfonts

fn=fn+1;figure(fn);clf;
    loglog(Nfc,L1fc./Lmfc,'k.','MarkerSize',20);
    hold on;grid on;
    loglog(Nfp,(0.873)*L1fp./Lmfp,'r.','MarkerSize',20);
    legend('$\lambda_{0}/\lambda_{N/2}$ (Fixed Charge)','$(1/0.873)\lambda_{0}/\lambda_{N/2}$ (Fixed Position)','Location','NorthWest');
    xlabel('N');
    %axis tight;
    set(gca,'ylim',[1.2,3.0]);
    setfonts

fn=fn+1;figure(fn);clf;
    semilogx(Nfc,Ffc_a,'k.','MarkerSize',20);
    hold on;grid on;
    semilogx(Nfc,Ffc_b,'k.','MarkerSize',5);
    semilogx(Nfp,Ffp_a,'r.','MarkerSize',20);
    semilogx(Nfp,Ffp_b,'r.','MarkerSize',5);
    legend('% $> 0.5$ (Fixed Charge)','% $\le 0.5$ (Fixed Charge)','% $> 0.5$ (Fixed Position)','% $\le 0.5$ (Fixed Position)','Location','NorthWest');
    xlabel('N');
    setfonts
    
fn=fn+1;figure(fn);clf;
    semilogx(Nfc,2*Afc_a,'k.','MarkerSize',20);
    hold on;grid on;
    semilogx(Nfc,2*Afc_b,'k.','MarkerSize',5);
    semilogx(Nfp,2*Afp_a,'r.','MarkerSize',20);
    semilogx(Nfp,2*Afp_b,'r.','MarkerSize',5);
    legend('$A > 0.5$ (Fixed Charge)','$A \le 0.5$ (Fixed Charge)','$A > 0.5$ (Fixed Position)','$A \le 0.5$ (Fixed Position)','Location','NorthWest');
    xlabel('N');
    setfonts

fn=fn+1;figure(fn);clf;
    semilogx(Nfc,xfc,'.','MarkerSize',20);
    grid on;hold on;
    title('$x_{i}$ when $\lambda_{i-1} > 0.5$ and $\lambda_{i} \le 0.5$');
    xlabel('N');
    setfonts
    
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
    legend('0.873$\lambda_{0} (Fixed Position)','$\lambda_{0} (Fixed Charge)','0.8$\lambda_{0} (Bonifirm Fit)','0.873$\lambda_{0} (Fixed Position Fit for N > 100)','Location','NorthWest');
    xlabel('N');
    setfonts
    
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
    legend('d$\lambda_{0}/dN (Fixed Position)','d$\lambda_{0}/dN (Fixed Charge)','Location','NorthWest');
    xlabel('N');
    setfonts
    