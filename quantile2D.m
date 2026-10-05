function [xdata,ydata,edata,n,qGebruik] = quantile2D(xdatain,ydatain,binstart,binwidth,binend)

binvect     = [binstart:binwidth:binend];
qUsable     = ~isnan(xdatain) & ~isnan(ydatain);
xdata       = quantile(xdatain(qUsable),binvect);

edges       = quantile(xdatain(qUsable),[binstart:binwidth:binend+binwidth]-binwidth/2);
[~,ind]     = histc(xdatain,edges);
qGebruik    = ind>0 & ind<length(edges);

ydata       = accumarray(ind(qGebruik),ydatain(qGebruik),[],@(x)median(x,'omitnan'));
ydata       = ydata';  % om dimensies hetzelfde te maken als xdata
[~,n]       = UniqueAndN(ind(qGebruik & ~isnan(ydatain)));

% MAD
edata       = accumarray(ind(qGebruik),ydatain(qGebruik),[],@(x)median(abs(x-median(x,'omitnan')),'omitnan'));
edata       = edata';  % om dimensies hetzelfde te maken als xdata

