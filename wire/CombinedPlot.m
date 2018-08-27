clear;

set(0,'DefaultFigureWindowStyle','docked')

xx = -0.99;

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
    
    Xfp{k}   = linspace(-1,1,Nfp(k));
    xdfp(k)  = Xfp{k}(find(Lfp{k}<=0.5,1));
    
    Nfp_a(k) = length(find(Lfp{k}>0.5));
    Nfp_b(k) = length(find(Lfp{k}<=0.5));
    Ffp_a(k) = Nfp_a(k)/Nfp(k); % Fraction above
    Ffp_b(k) = Nfp_b(k)/Nfp(k); % Fraction below
    Afp_a(k) = (1/2)*sum(Lfp{k}(Lfp{k}>0.5))/Nfp_a(k); % Area above
    Afp_b(k) = (1/2)*sum(Lfp{k}(Lfp{k}<=0.5))/Nfp_b(k); % Area above

    I = find(Xfp{k}>xx);
    %FFfp(k) = Lfp{k}(I(1));
    Ix = -1+ceil(sqrt(Nfp(k)));
    FFfp(k) = mean(Lambda{k}(Ix))*Nfp(k)/2;

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

        ds     = d*(Nfc(k));
        xd     = ds/2+Xfc{k}(1:end-1); % Center point between each charge.
        dave(k) = 2.0;
        xdfc(k) = Xfc{k}(find(ds>dave(k),1)); % First location where delta <= dave.
        Nfc_a(k) = length(find(ds>dave(k)));
        Nfc_b(k) = length(find(ds<=dave(k)));
        Ffc_a(k) = Nfc_a(k)/Nfc(k); % Fraction above
        Ffc_b(k) = Nfc_b(k)/Nfc(k); % Fraction below
        Afc_a(k) = (1/2)*sum(ds(ds>dave(k)))/Nfc_a(k); % Area above
        Afc_b(k) = (1/2)*sum(ds(ds<=dave(k)))/Nfc_b(k); % Area above

        %I = find(Xfc{k}<=-0.99);
        I = find(Xfc{k}>xx);
        FFfc(k) = 1./(Nfc(k)*mean( d(1:-1+ceil(0.5*sqrt(Nfc(k))) )));
        
        % Assume total charge is 1.
        % Lamda = dq/dx = (1/N)/(x_{i+1}-x_{i})
        Lfc{k}   = 1./(Nfc(k)*d);
        if (0)
            Lave(k)  = mean(Lfc{k});
            Lc       = 0.5;
            xfc(k)   = X{i}(find(Lfc{k}<=Lc,1)); % First location where Lambda <= 0.5.
            Nfc_a(k) = length(find(Lfc{k}>Lc));
            Nfc_b(k) = length(find(Lfc{k}<=Lc));
            Ffc_a(k) = Nfc_a(k)/Nfc(k); % Fraction above
            Ffc_b(k) = Nfc_b(k)/Nfc(k); % Fraction below
            Afc_a(k) = (1/2)*sum(Lfc{k}(Lfc{k}>Lc))/Nfc_a(k); % Area above
            Afc_b(k) = (1/2)*sum(Lfc{k}(Lfc{k}<=Lc))/Nfc_b(k); % Area above
        end
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
    set(gca,'xlim',[-1.02,0.02])
    ylabel('$\lambda$ (Fixed Charge)');
    %z = 0.01;
    %x = [-1+z:z:1-z];
    %y = 1./sqrt(1-(x).^2);
    %plot(x,y-min(y)+0.45,'k','LineWidth',2);
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
        d = Lfc{i};
        dx = diff(Xfc{i});
        x = dx/2+Xfc{i}(1:end-1); % Center point between each charge.
        plot(x,d,'Marker','.','LineStyle','-');
    end
    legend off
    title(Is);
    xlabel('x');
    set(gca,'xlim',[-0.95,-0.55])
    ylabel('$\lambda$ (Fixed Charge)');
    setfonts    
    
fn=fn+1;figure(fn);clf;grid on;hold on;
    for i = I
        d = diff(Xfc{i});
        x = d(1)/2+Xfc{i}(1:end-1); % Center point between each charge.
        plot(x,d*length(Xfc{i}));
    end
    legend off
    title(Is);
    xlabel('x');
    set(gca,'xlim',[-1.02,0.02])
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
    semilogx(Nfc,xdfc,'.','MarkerSize',20);
    grid on;hold on;
    semilogx(Nfp,xdfp,'.','MarkerSize',20);
    legend('$x_{i}$ when $\Delta_{i-1} > 2$ and $\Delta_{i} \le 2$ (Fixed Charge)','$x_{i}$ when $\lambda_{i-1} > 0.5$ and $\lambda_{i} \le 0.5$ (Fixed Position)');
    ylabel('$x_c$');
    set(gca,'ylim',[-1,-0.6])
    xlabel('N');
    setfonts
    
fn=fn+1;figure(fn);clf;
    loglog(Nfc,xdfc,'.','MarkerSize',20);
    grid on;hold on;
    loglog(Nfp,xdfp,'.','MarkerSize',20);
    legend('$x_{i}$ when $\Delta_{i-1} > 2$ and $\Delta_{i} \le 2$ (Fixed Charge)','$x_{i}$ when $\lambda_{i-1} > 0.5$ and $\lambda_{i} \le 0.5$ (Fixed Position)');
    ylabel('$x_c$');
    set(gca,'ylim',[-1,-0.6])
    xlabel('N');
    setfonts
    
if (0)    
    p = polyfit(log10(Nfc(1000:end-1)),log10(-diff(xdfc(1000:end))),1)
    p = polyfit(log10(Nfc(100:1000-1)),log10(-diff(xdfc(100:1000))),1)

    y = polyval(p,log10(Nfc(100:end)));
    fn=fn+1;figure(fn);clf;
        loglog(Nfc(100:end),-10.^y,'b','LineWidth',5);
        grid on;hold on;
        loglog(Nfc(1:end-1),diff(xdfc),'g.','MarkerSize',3);
        legend('Fit for N=100+','diff($\Delta$)','Location','NorthWest');
        xlabel('N');
        setfonts
end
    
fn=fn+1;figure(fn);clf;
    semilogx(Nfp,Afp_a,'r.','MarkerSize',20);
    hold on;grid on;
    semilogx(Nfp,Afp_b,'r.','MarkerSize',5);
    legend('$A > x_c$ (Fixed Position)','$A \le x_c$ (Fixed Position)','Location','NorthWest');
    xlabel('N');
    setfonts
    
fn=fn+1;figure(fn);clf;
    semilogx(Nfc,Afc_a,'k.','MarkerSize',20);
    hold on;grid on;
    semilogx(Nfc,Afc_b,'k.','MarkerSize',5);
    legend('$A > x_c$ (Fixed Charge)','$A \le x_c$ (Fixed Charge)','Location','NorthEast');
    xlabel('N');
    setfonts
 
fn=fn+1;figure(fn);clf;
    semilogx(Nfc,Ffc_a,'k.','MarkerSize',20);
    hold on;grid on;
    semilogx(Nfc,Ffc_b,'k.','MarkerSize',5);
    semilogx(Nfp,Ffp_a,'r.','MarkerSize',20);
    semilogx(Nfp,Ffp_b,'r.','MarkerSize',5);
    legend('Fraction $> x_c$ (Fixed Charge)','Fraction $\le x_c$ (Fixed Charge)','Fraction $> x_c$ (Fixed Position)','Fraction $\le x_c$ (Fixed Position)','Location','NorthWest');
    xlabel('N');
    setfonts

fn=fn+1;figure(fn);clf;
    semilogx(Nfc,FFfc,'.','MarkerSize',10);
    grid on;hold on;
    semilogx(Nfp,FFfp,'.','MarkerSize',15)
    set(gca,'ylim',[0,0.14])
    legend('Fixed Charge','Fixed Position');
    title(sprintf('Fraction charges with x < %.2f',xx));
    xlabel('N');

fn=fn+1;figure(fn);clf;
    loglog(Nfc,FFfc,'.','MarkerSize',10);
    grid on;hold on;
    loglog(Nfp,FFfp,'.','MarkerSize',15)
    %set(gca,'ylim',[0,0.14])
    legend('Fixed Charge','Fixed Position');
    %ylabel('N\Delta_0')
    title(sprintf('$\\lambda\\simeq$ %.2f',xx));
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
    