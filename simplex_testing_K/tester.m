function [] = tester
xc_int = [0;0;0];  
xc = randn(2,1);    xc = xc./norm(xc,2)*(1E-1+1E-1*rand(1));    
xc = [xc(1)^2+xc(2)^2;xc(1);xc(2)];
xc_front = randn(2,2);  
xc_front(:,1) = xc_front(:,1)./norm(xc_front(:,1),2)*(1E-1+1E-1*rand(1));
xc_front(:,2) = xc_front(:,2)./norm(xc_front(:,2),2)*(1E-1+1E-1*rand(1));
xc_front = [zeros(1,2);xc_front] + [0;xc(2:3)];
xc_front(1,1) = xc_front(2,1)^2 + xc_front(3,1)^2;
xc_front(1,2) = xc_front(2,2)^2 + xc_front(3,2)^2;
xc_h = min([norm(xc,2),norm(xc-xc_front(:,1),2),norm(xc-xc_front(:,2),2)]);
[x_new,R2frame,yframe] = grow_simplex(xc,xc_front,xc_int,xc_h);
if isempty(x_new)
    tester
else
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot(2,2,1)
    plot_para_frame(x_new,xc);
    plotvect(xc(2:3),[0;0],'red')
    xlabel('a','FontSize',20)
    ylabel('b','FontSize',20)
    title('Parameter frame')
    subplot(2,2,2)
    plot_R2_frame(R2frame)
    title('Normalized tangent space projection on R2')
    subplot(2,2,3)
    plot_tangent_frame(yframe);
    title('Normalized tangent space - not projected.')
    subplot(2,2,4)
    [~,m] = size(x_new);
    % Rotate coordinates so that (a,b) are in the plane.
    x_new = [x_new(2:end,:);x_new(1,:)];
    plot3([0,x_new(1,1)],[0,x_new(2,1)],[0,x_new(3,1)],'red')
    hold on
    for k=1:m-2
        h = patch('Faces',1:3,'Vertices',[x_new(:,1)';x_new(:,k+1)';x_new(:,k+2)']);
        set(h,'FaceColor','cyan','EdgeColor','b','LineWidth',2,'FaceAlpha',0.5);
        xlabel('a','FontSize',20)
        ylabel('b','FontSize',20)
        zlabel('x','FontSize',20)
        title('Manifold')
    end
end
view(3)
end