function [] = plot_tangent_frame(x)
[~,m]=size(x);
plotvect3(x(:,1),[0;0;0],'k');
hold on
plotvect3(x(:,2),[0;0;0],'r');
plotvect3(x(:,3),[0;0;0],'k');
for j=4:m
    plotvect3(x(:,j),[0;0;0],'blue')
end
end