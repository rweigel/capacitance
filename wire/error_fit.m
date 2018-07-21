function [err,lfit] = error_fit(p,n,ldata)

%lfit = log(n).*(p(1) + p(2)./log(n) + p(3)./(log(n)).^2);

I = find(n>100);
n = n(I);
ldata = ldata(I);
lfit = p(1) + p(2).*log(n);

err = sum((lfit - ldata).^2);