
            if (0)
                    % Calculate 2-parameter model for n particles.
                    options = optimset('TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',1000);
                    p1 = fminsearch(@(p) pushcode_err(p,n),p1,options);
                    [Fss,x] = pushcode_err(p1,n);
                    Fss
                end

                
                % Calculate 2-parameter model for n-1 particles.
                options = optimset('TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',1000);
                p = fminsearch(@(p) pushcode_err(p,n-1),p,options);
                [Fss,xz] = pushcode_err(p,n-1);
                %Fss
                
                % Difference between last calculated and 2-parameter model
                % for n-1 particles ("model error").
                lx = X{n-1}-xz;

                % Interpolate model error on to n point grid.
                lxi = interp1(linspace(-1,1,n-1),lx,linspace(-1,1,n),'pchip');
                
                % Calculate 3-parameter model for n particles
                % Model is 
                %   2-parameter model + const3*(interpolated model error)
                options = optimset('TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',1000);
                p = fminsearch(@(p) pushcode_err(p,n,lxi),p,options);
                [Fss,x] = pushcode_err(p,n,lxi);
     
clear;
set(0,'DefaultFigureWindowStyle','docked');
figure(1);clf;
figure(2);clf;

if true && exist('pushcode.mat','file')
    load pushcode;
    no = length(X{end});
    x  = X{no};
else
    no = 3; % Initial # of particles
    x  = linspace(-1,1,no);
end

dn = 1; % Step in number of particles
nf = 10000; % Final number of particles
po = [9.75,1000,-1e-3];
for n = no:dn:nf
    tic;
    alpha  = 1;
    t      = 0;
    t_stop = 500;
    Fss_stop = 1e-5;%n*eps*1e10;
    Fss      = Fss_stop*10;
    
    if exist('X','var')
        if (n > 10)
            if ~exist('p','var')
                p = [0.91,5,1];
            end
            %options = optimset('TolFun',1e-9,'TolX',1e-9,'Display','iter','MaxFunEvals',50000);
            %p = fminsearch(@(p) pushcode_err(p,n),p,options);
            %[Fmse,x] = pushcode_err(p,n);
            if (1)                 
                options = optimset('TolFun',1e-9,'TolX',1e-9,'MaxFunEvals',1000);
                p = fminsearch(@(p) pushcode_err(p,n-1),po,options);
                [Fss,xz] = pushcode_err(p,n-1);

                lx = X{n-1}-xz;
                lxi = interp1(linspace(-1,1,n-1),lx,linspace(-1,1,n),'pchip');
                p = fminsearch(@(p) pushcode_err(p,n,lxi),p,options);
                [Fss,x] = pushcode_err(p,n,lxi);
                xo = x;
            end
            
            if (0)
                xo = linspace(-1,1,n-1);
                D = diff(x);
                p = polyfit(linspace(-1,1,length(D)),D,9);
                tmpx = polyval(p,xo);
                off = [0,cumsum(tmpx)];
                xo = linspace(-1,1,n);
                x = -1+2*off/off(end);
            end
            
        else
            xo = linspace(-1,1,n);            
            xolast = linspace(-1,1,length(x));
            x = interp1(xolast,x-xolast,xo,'spline')+xo;
            %x = interpft(x-xolast,n)+xo;
            %x = interp1(xolast,x-xolast,xo)+xo;
            %x = interp1(xolast,x,xo);
        end
    else
        x  = linspace(-1,1,n);
        xo = x; % Starting guess positions
    end
    
    N = length(x);
    I = [1:N];
    
    while Fss > Fss_stop && t < t_stop 
        t = t+1;
        F = zeros(1,N);
        for i = 1:N
            for j = 1:N
                if i ~= j
                    F(i) = F(i) + sign(x(i)-x(j))/(x(i)-x(j))^2;
                end
            end
        end
        Fss = sqrt(sum(F(2:end-1).^2))/(length(F)-2); % External force is not zero.
        %fprintf('t = %03d; a = %.2e; Fss = %.2e\n',t,alpha,Fss);
        if (Fss < Fss_stop)
            fprintf('n = %02d; t = %03d; alpha = %.2e; Fss = %.2e Fstop = %.2e\n',n,t,alpha,Fss,Fss_stop);
            break;
        end

        r = 1;
        xt = x;
        if (r == 1)
            It = 0*I;
        end
        cont = 1;
        while cont
            found = 0;
            for i = 1:N
                xt(i) = x(i) + alpha*F(i);
                if xt(i) < -1
                    xt(i) = -1; 
                end
                if xt(i) > 1
                    xt(i) = 1;
                end
                if (0)
                    if (i>2)
                        if ~found && (~(xt(i-2) < xt(i-1) && xt(i-1) < xt(i)) || any(abs(diff(xt))<10*eps))
                            %fprintf('x\n')
                            found = 1;
                            alpha = alpha/2;
                            r = r+1;
                            break
                        end
                    end
                end
            end
            if (0)
                if (found == 0)
                    break
                end
            end
            if (1)
                [jnk,It] = sort(xt);
                cont = Ferr(xt) > Fss || any(It - I) || any(abs(diff(xt))<10*eps);
                if cont
                    %fprintf('y\n')
                    alpha = alpha/2;
                    r = r+1;
                end
            end
        end
        x = xt;        
        if mod(t,100) == 0 || t == 1 || t == t_stop
            fprintf('n = %02d; t = %03d; alpha = %.2e; Fss = %.2e Fstop = %.2e\n',n,t,alpha,Fss,Fss_stop);
            if (t==1);
                figure(1);clf;hold on;grid on;
            else
                set(0,'CurrentFigure',1);
            end
            plot(xo,(x-xo)/max(x),'b.');hold on;
            %plot(xo(2:end-1)',(x(2:end-1)'-xo(2:end-1)')./xo(2:end-1)','r.');
            drawnow;
            if (n==no);
                figure(2);
            else
                set(0,'CurrentFigure',2);
            end
            plot(x,x-xo,'-');
            drawnow;
        end
       
    end
    set(0,'CurrentFigure',2);clf;hold on;grid on;
    if exist('p','var')
        P{n} = p;
    end
    X{n} = x;
    Xo{n} = xo;
    Y{n} = x-xo; % Equilib. pos. - linear pos.
    T{n} = toc;
    A{n} = alpha;

    save pushcode X Xo T A
      
end