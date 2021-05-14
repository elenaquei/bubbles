function [] = plot_para_frame(x,xc)
x = x(2:end,:);
[~,m]=size(x);
plotvect(x(:,1),xc(2:3),'*r');
hold on
plotvect(x(:,2),xc(2:3),'k');
plotvect(x(:,m),xc(2:3),'k');
for j=3:m-1
    plotvect(x(:,j),xc(2:3),'blue')
end
axis equal
end

