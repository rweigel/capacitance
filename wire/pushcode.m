clear;
set(0,'DefaultFigureWindowStyle','docked');
figure(1);clf;
figure(2);clf;
mkdir('tmp');

if false && exist('pushcodec.mat','file')
    load pushcodec.mat;
    no = length(X{end});
    % Re-do last iteration for comparison
    for i = length(X)-1:-1:1
        if ~isempty(X{i})
            break;
        end
    end
    no = i;
    x  = X{no};
else
    no = 3; % Initial # of particles
    x  = linspace(-1,1,no);
end
        
dn = 1; % Step in number of particles
nf = 5000; % Final number of particles
nf = 100;

system('make');

DataM = struct();
DataC = struct();

for n = no:dn:nf
    alpha  = 1;
    t      = 0;
    t_stop = 5000;
    Fss_stop = 1e-5;%n*eps*1e10;
    Fss      = Fss_stop*10;
    tic;
    if exist('X','var')
        if (n > 30)
            if ~exist('p','var')
                po = [9.75,1000];
                p  = [9.75,1000,0.001];
                %p = [0.91,5,1];
                %p = [0.91,5];
            end
            if (0)
                % (Worse performance)
                % Calculate 2-parameter model for n particles.
                options = optimset('TolFun',1e-9,'TolX',1e-9,'Display','iter','MaxFunEvals',50000);
                p = fminsearch(@(p) pushcode_err(p,n),p,options);
                [Fmse,x] = pushcode_err(p,n);
                Fmse
            end
            if (1)
                options = optimset('TolFun',1e-2,'TolX',1e-2,'MaxFunEvals',1000);
                % Calculate 2-parameter model for n-1 particles.
                po = fminsearch(@(p) pushcode_err(p,n-dn),po,options);
                [Fss,xz] = pushcode_err(po,n-dn);
                %Fss
                %p'

                % Difference between last calculated and 2-parameter model
                % for n-1 particles.
                lx = X{n-dn}-xz;

                % Interpolate model error on to n point grid.
                lxi = interp1(linspace(-1,1,n-dn),lx,linspace(-1,1,n),'pchip');
                p = fminsearch(@(p) pushcode_err(p,n,lxi),p,options);
                [Fss,x] = pushcode_err(p,n,lxi);
            end
            xo = x;
        else
            xo = linspace(-1,1,n);            
            xolast = linspace(-1,1,length(x));
            x = interp1(xolast,x-xolast,xo,'spline')+xo;
            % Worse performance:
            %x = interpft(x-xolast,n)+xo;
            %x = interp1(xolast,x-xolast,xo)+xo;
            %x = interp1(xolast,x,xo);
        end
    else
        x  = linspace(-1,1,n);
        xo = x; % Starting guess positions
    end
    tfit = toc;
    
    N = length(x);
    I = [1:N];
    
    fid = fopen('./tmp/pushcode_xo.txt','w');
    fprintf(fid,'%.16f\n',x');
    fclose(fid);
    tic;
    system(sprintf('./pushcode %d %d',n,t_stop));
    Z = load('./tmp/pushcode_xf.txt','r');
    
    DataC.Time(n)  = toc;
    DataC.TimeFit(n) = tfit;
    DataC.Niter(n) = Z(1);
    DataC.FSS(n)   = Z(2);
    DataC.Alpha(n) = Z(3);
    DataC.X{n}     = Z(4:end)';
    DataC.Xo{n}    = xo;

    fprintf('c: n = %03d; T = %03d; alpha = %.2e; Fss = %.2e; tcode = %.2e; tfit = %.2e\n',n,DataC.Niter(n),DataC.Alpha(n),DataC.FSS(n),DataC.Time(n),tfit);
    save ./mat/pushcodec DataC
    figure(1);clf;hold on;grid on;
        plot(xo,(DataC.X{n}-xo),'b.');
        xlabel('x');
        legend('x_f-x_o','Location','NorthWest');
        drawnow;

    if ~(n < 50 || n == 100 || n == 1000)
        continue;
    end
    
    tic;
    while Fss > Fss_stop && t < t_stop
        alpha = 1;
        t = t+1;
        [Fss,F] = Ferr(x);
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
        Fss = Ferr(xt);               
    end
    
    DataM.Time(n)    = toc;
    DataM.TimeFit(n) = tfit;
    DataM.Niter(n) = t;
    DataM.FSS(n)   = Fss;
    DataM.Alpha(n) = alpha;
    DataM.X{n}     = x;
    DataM.Xo{n}    = xo;

    fprintf('m: n = %03d; T = %03d; alpha = %.2e; Fss = %.2e; tcode = %.2e; tfit = %.2e\n',n,DataM.Niter(n),DataM.Alpha(n),DataM.FSS(n),DataM.Time(n),tfit);
    save ./mat/pushcodem.mat DataM

    figure(2);clf;hold on;grid on;
        plot(xo,DataM.X{n}-xo,'b.');
        xlabel('x');
        legend('x_f-x_o','Location','NorthWest');
        drawnow;
    
end