clear
% Force method

% For n odd, get ill-conditioned matrix For n even, get different results
% depending on alpha. Get results that match Griffiths when alpha = 1 and N
% = 4 and 6 (his n=2 and 4). Plots for N=10, ... don't match his plots
% exactly, however (my end values are higher).

alpha = 0.95;
k = 1;
%for N = 10:100:10000;
No = 2.^[3:18];
for i = 1:length(No);
    N = No(i);    
    if (0) % Set spacing to decrease near ends.
        if mod(N,2) == 0
            x(1) = alpha/2;
            Nr = N/2-1;
        else
            x(1) = 0;
            Nr = floor(N/2);
        end
        for i=1:Nr
            x(i+1)=x(i)+(alpha)^i;
        end
        x = x/(2*max(x));
        if mod(N,2) == 0
            x = [-fliplr(x(1:end)),x];
        else
            x = [-fliplr(x(2:end)),x];
        end
    else
        % Linear spacing
        x = linspace(-0.5,0.5,N);
    end

    Z = zeros(N);
    tic
    for i = 1:N
        if i > 13 && mod(i,100) == 0
           fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b',i,N);
           fprintf('%d/%d',i,N);
        end
        for j = 1:N
            if (i ~= j)
                d = x(i)-x(j);
                Z(i,j) = sign(d)/d^2;
            end
        end
    end
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b',i,N);

    b = [zeros(size(Z,1),1);1]; % Total charge
    Z(end+1,:) = ones(1,size(Z,2));
    Z(1,end+1) = 1; % Wall force
    Z(end-1,end) = -1; % Wall force
    t0 = toc;
    
    tic
    l = Z\b;
    t1 = toc;
    F1 = l(end);
    l1 = l(1:end-1);

    L{k} = l1;
    R1(k) = l1(1)/l1(end/2);
    save FixedPosition.mat L R1
    
    if (0)
    tic
    [u,s,v] = svd(Z);
    l2 = v(1:end-1,end);
    R2(k) = l2(1)/l2(end/2);

    F2 = l2(end,end);
    t2 = toc;
    end
    fprintf('N = %d; R1 = %.2f; t0 = %.1e; t1 = %.1e;\n',N,R1(k),t0,t1);
    %fprintf('N = %d; R1 = %.2f; R2 = %.2f; t1 = %.1e; t2 = %.1e\n',N,R1(k),R2(k),t1,t2);
    k = k+1;
end


figure(1);clf;grid on;hold on;
    %plot(x,l1,'LineWidth',2);
    %plot(x,l2*N/2,'LineWidth',2);
    plot(x,l1*N/2,'.','MarkerSize',20);
    %plot(x,l2*N/2,'.','MarkerSize',20);
    plot([-0.5,0.5],[0,0],'k')
    xlabel('x''');
    %ylabel('$$\frac{kq}{LV_o}$$','Interpreter','Latex','Rotation',0,'HorizontalAlignment','Right')
    box on;
    legend(sprintf('\\lambda''_1; q''_1=%.2g',sum(l1)));%,...
        %sprintf('\\lambda''_2; q''_1=%.2g',sum(l2)));
saveplots('./figures/Capacitance_Force_Method_Density');
