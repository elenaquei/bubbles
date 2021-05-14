function [] = plot_R2_frame(x)
[~,m]=size(x);
plotvect(x(:,1),[0;0],'k');
hold on
plotvect(x(:,2),[0;0],'r');
plotvect(x(:,3),[0;0],'k');
for j=4:m
    plotvect(x(:,j),[0;0],'blue')
end
axis equal
end

