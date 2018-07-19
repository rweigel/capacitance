function [Fss,F] = Ferr(x)

F = zeros(1,length(x));
for i = 1:length(x)
    F(i) = 0;
    for j = 1:length(x)
        if i ~= j
            F(i) = F(i) + sign(x(i)-x(j))/(x(i)-x(j))^2;
        end
    end
end

% Boundary forces are not zero.
Fss = sqrt(sum(F(2:end-1).^2))/(length(F)-2);