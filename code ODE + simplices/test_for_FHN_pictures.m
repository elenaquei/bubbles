% testing for pretty pictures in FHN

% Crude experiments with MATLAB integrator: It looks like there is a singular
% Hopf bifurcation near (?, I, ?, ?) = (0.1, 0.4, 0.305, 1).

epsilon = linspace( 0.2,  0.33 , 50);
I = linspace(0.125, 0.5 , 50);
figure(1)
figure(2)
amplitude = 0*I;

init_coord  = [0.4491    0.4829]; 
approx_period = 10*14.7;
        
for j = 1:length(epsilon)
    epsilon_fix = epsilon(j);
    for i = 1:length(I)
        alpha = 0.1;
        I_fix = I(i);
        % epsilon_fix = 0.2480;
        % testing_epsilon_fix = .180;
        % gamma = 1; % will not be coded for simplicity
        
        % right hand side
        f=@(t,x)[x(1)^2 - alpha * x(1) - x(1)^3 + alpha * x(1)^2 - x(2) + I_fix ;
            epsilon_fix * x(1) - epsilon_fix * x(2)];
        
        % forward integration
        [~, yout] = ode23s(f,linspace(0,10*approx_period,200),init_coord);
        
        init_coord2= yout(end,:);
        
        [~, yout] = ode23s(f,linspace(0,approx_period,200),init_coord2);
        
        %plot(yout(:,1),yout(:,2))
        
        amplitude(i) = max(yout(:,1)) - min(yout(:,1));
        orbit(i,:) = [max(yout(:,1)) , min(yout(:,1))];
    end
    
    figure(1)
    plot(I, amplitude, 'LineWidth', 2);
    axis([I(1) I(end) 0 0.8])
    %axis([I(1) I(end) -0.2 0.8])
    set(gca,'FontSize',20)
    s = ['$\epsilon$ = ',num2str(epsilon_fix)];
    title(s,'Interpreter','Latex', 'FontSize', 24);
    xlabel('$I$','Interpreter','Latex', 'FontSize', 24);
    ylabel('amplitude','Interpreter','Latex', 'FontSize', 24);
    drawnow
    F(j) = getframe(gcf) ;
    drawnow
    figure(2)
    plot(I, orbit(:,1), 'LineWidth', 2);
    hold on
    plot(I, orbit(:,2), 'LineWidth', 2);
    hold off
    axis([I(1) I(end) -0.2 0.8])
    set(gca,'FontSize',20)
    s = ['$\epsilon$ = ',num2str(epsilon_fix)];
    title(s,'Interpreter','Latex', 'FontSize', 24);
    xlabel('$I$','Interpreter','Latex', 'FontSize', 24);
    ylabel('minimum and maximum','Interpreter','Latex', 'FontSize', 24);
    drawnow
    G(j) = getframe(gcf) ;
    drawnow
end

% error('stop here for now')

% create the video writer with 1 fps
writerObj = VideoWriter('amplitude_vs_I_and_epsilon.mp4','MPEG-4');
writerObj.FrameRate = 5;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(F)
    % convert the image to a frame
    frame = F(i) ;
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);

% create the video writer with 1 fps
writerObj = VideoWriter('minimum_and_maximum_vs_I_and_epsilon.mp4','MPEG-4');
writerObj.FrameRate = 5;
% set the seconds per image
% open the video writer
open(writerObj);
% write the frames to the video
for i=1:length(G)
    % convert the image to a frame
    frame = G(i) ;
    writeVideo(writerObj, frame);
end
% close the writer object
close(writerObj);
