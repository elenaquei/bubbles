function surface_plot_lorenz84(iter1, iter2, iter3)
% function surface_plot_lorenz84(iter1, iter2, iter3)
%
% this function assumes the main lorenz84_cont_2 has been run successfully
% and builds on the files that that main stored to plot the orbits of the
% lorenz 84 model between two Hopf bifurcations


if nargin == 0 || ~exist('iter1','var')
    iter1=60;
end
if nargin < 2 || ~exist('iter2','var')
    iter2=20;
end
if nargin<3 || ~exist('iter3','var')
    iter3=33;
end

number_of_Hopf1 = 1:iter1; % 10
number_of_not_Hopf = iter1 + (1:iter2); % 45
number_of_Hopf2 = iter1 + iter2 + (1:iter3); %19
names_savefiles = cell(max([number_of_Hopf1,number_of_Hopf2]),1);
names_savefiles_middle = cell(numel(number_of_not_Hopf),1);

max_n_nodes = 0;


for i = number_of_Hopf1
    names_savefiles{i} = strcat('lorenz84_close_to_Hopf',sprintf('%i',i));
    load(names_savefiles{i})
    max_n_nodes = max(max_n_nodes, x_n.nodes);
end
iter = 1;
for i = number_of_not_Hopf
    names_savefiles_middle{i} = strcat('Hopf_lorenz84_after_Hopf',sprintf('%i',iter));
    iter = iter+1;
    load(names_savefiles_middle{i})
    max_n_nodes = max(max_n_nodes, x_n.nodes);
end
iter = 1;
for i = number_of_Hopf2
    names_savefiles{i} = strcat('Hopf_lorenz84_to_second_Hopf',sprintf('%i',iter));
    load(names_savefiles{i})
    max_n_nodes = max(max_n_nodes, x_n.nodes);
    iter = iter+1;
end

all_Fouriers = zeros(iter1+iter2+iter3,4, 2*max_n_nodes+1);
lambda = zeros(iter1+iter2+iter3,1);

for i = number_of_Hopf1
    load(names_savefiles{i})
    x_n = reshape_Xi(x_n,max_n_nodes);
    x_normal = rescale_Hopf(x_n);
    all_Fouriers(i,:,:) = x_normal.vector;
    lambda(i) = x_n.scalar(2);
    %disp(x_n.scalar(2))
end

for i = number_of_not_Hopf
    load(names_savefiles_middle{i})
    x_n = reshape_Xi(x_n,max_n_nodes);
    all_Fouriers(i,:,:) = x_n.vector;
    lambda(i) = x_n.scalar(2);
    %disp(x_n.scalar(2))
end

for i = number_of_Hopf2
    clear saddle_x0_stored
    load(names_savefiles{i})
    
    if exist('saddle_x0_stored', 'var')
        last_element = i;
        number_of_Hopf2 = iter1+iter2 : last_element;
        
        x_n_possible1 = saddle_x0_stored{1};
        x_n_possible2 = saddle_x1_stored{1};
        
        if x_n_possible1.scalar(2)>x_n_possible2.scalar(2)
            x_n  = x_n_possible2;
        else
            x_n  = x_n_possible1;
        end
    end
    
    x_n = reshape_Xi(x_n,max_n_nodes);
    x_normal = rescale_Hopf(x_n);
    all_Fouriers(i,:,:) = x_normal.vector;
    lambda(i) = x_n.scalar(2);
    
    %disp(x_n.scalar(2))
    if exist('saddle_x0_stored', 'var')
        break
    end
end
all_Fouriers = all_Fouriers(1:last_element,:,:);
lambda = lambda(1:last_element);

Dims = [1,3;
    2,4];

% view1 = [-37.5, -80.7];
% view2 = [30,20.5];

view1=[-20,-20];
view2=[50,5];
volume{1}=[-0.0175,0.0175,0.85,1.3,-0.15,0.75];
volume{2}=[-0.0175,0.0175,-0.39,0.39,-0.5,0.5];

%color = colormap(parula);
%index_step = floor(size(color,1)/max(number_of_Hopf2));
%colors = color(1:index_step:size(color,1),:);


clf;

for dim_loc = [1,2]
    subplot(1,2,dim_loc)
    dims = Dims(dim_loc,:);
    
    load(names_savefiles{1})
    %plot_Hopf(x_0,dims,'k.','MarkerSize',40)
    plot_numerical_equilibria(x_0.scalar(2), x_0.scalar(4:end), dims);
    hold on
    
     
    plot_surface(all_Fouriers(:,dims,:),lambda, iter1, iter2)
    
   view(view1(dim_loc), view2(dim_loc)); % set point of view for better sight

   load(names_savefiles{number_of_Hopf2(end)})
    x_Hopf_0 = saddle_x0_stored{1};
    x_Hopf_1 = saddle_x1_stored{1};
    amplitudes = [x_Hopf_0.scalar(3), x_Hopf_1.scalar(3)];
    [~,index] = min(amplitudes);
    % find the linear interpolation with amplitude = 0
    if index == 1
        x_Hopf = x_Hopf_0;
    else
        x_Hopf = x_Hopf_1;
    end
    %plot_Hopf(x_Hopf,dims,'r.','MarkerSize',40)
    plot_numerical_equilibria(x_Hopf.scalar(2), x_Hopf.scalar(4:end), dims);
    %plot_Hopf(x_Hopf,dims,'*','color',colors(end,:),'MarkerSize',10)
    xlabel('$\mu$','interpreter','Latex');
    ylabel(sprintf('$u_%i$',dims(1)),'Interpreter','Latex');
    zlabel(sprintf('$u_%i$',dims(2)),'Interpreter','Latex');
%    ylabel(sprintf('dimension %i',dims(1)),'interpreter','tex','Fontsize',15)
%    zlabel(sprintf('dimension %i',dims(2)),'interpreter','tex','Fontsize',15)
    axis(volume{dim_loc});
    ax = gca;
    ax.FontSize = 15;
   
end

set(gcf,'position',[100,200,1000,400])
hold off
% set(gcf,'PaperUnits','centimeters')
% set(gcf,'PaperSize',[31 15])
%printdir='';
%print([printdir,'balloons'],'-dpdf','-painters')

%printdir='../../tex/';
%print([printdir,'balloons'],'-djpeg','-painters','-r600')

end


function x_Xivec_rescaled = rescale_Hopf(x_Xivec)
x_Xivec_rescaled = x_Xivec;
amplitude = x_Xivec.scalar(3);
x_Xivec_rescaled.vector = (amplitude) * x_Xivec_rescaled.vector;
x_Xivec_rescaled.vector(:,x_Xivec.nodes+1) = x_Xivec.scalar(4:end)'+...
    x_Xivec_rescaled.vector(:,x_Xivec.nodes+1);
end




function [] = plot_surface(fouriers, parameters, iter1, iter2)

points=size(fouriers,3);
nt=100;
time=linspace(0,2*pi,nt);
nodes = (points-1)/2;

% number of solutions
N=size(fouriers,1); 


h=linspace(-1,1,N);

x=zeros([N,nt]);
y=zeros([N,nt]);
z=zeros([N,nt]);

for k=1:N
    % defined the N periodic solutions; 
    x(k,:)=parameters(k);
    for i = 1:length(time)
        y(k,i)=sum(squeeze(fouriers(k,1,:)).'.*exp(1i*[-nodes:nodes]*time(i)));
        z(k,i)=sum(squeeze(fouriers(k,2,:)).'.*exp(1i*[-nodes:nodes]*time(i)));
    end
    %y(k,:)=r(k)*cos(time);
    %z(k,:)=r(k)*sin(time);
end

if any(any(imag(y)>10^-10)) || any(any(imag(z)>10^-10))
    error('There is somethign fishy with the solutions')
end
y = real(y);
z = real(z);

h=surf(x,real(y),real(z));
set(h,'FaceAlpha',0.6)
shading interp
hold on
%view(-20,45);
h.CData = x;
colormap winter

for k=floor(linspace(1,N,10))
    if k<= iter1
        plot3(x(k,:),y(k,:),z(k,:),'k');
    elseif k<=iter1+iter2
        plot3(x(k,:),y(k,:),z(k,:),'k');
    else
        plot3(x(k,:),y(k,:),z(k,:),'k');
    end
end

%plot3(x(1,1),y(1,1),z(1,1),'k.','MarkerSize',30);
%plot3(x(N,1),y(N,1),z(N,1),'r.','MarkerSize',30);

end

function plot_Hopf(x_Xivec, dims,varargin)
x_Xivec_rescaled = x_Xivec;
amplitude = x_Xivec.scalar(3);
x_Xivec_rescaled.vector = (amplitude) * x_Xivec_rescaled.vector;
x_Xivec_rescaled.vector(:,x_Xivec.nodes+1) = x_Xivec.scalar(4:end)'+...
    x_Xivec_rescaled.vector(:,x_Xivec.nodes+1);
y = time_series_loc(x_Xivec_rescaled);
y = real(y);
lambda = x_Xivec.scalar(2).*ones(size(y,1),1);
plot3(lambda,mean(y(:,dims(1))),mean(y(:,dims(2))),varargin{:})


end


function [xt1,t1] = time_series_loc(X,X2)
if abs(X.scalar(1))>0.1
    npoints=200*ceil(1/real(X.scalar(1)));
    if nargin == 2
        npoints=max(200*ceil(1/real(X.scalar(1))),200*ceil(1/real(X2.scalar(1))));
    end
    t1=linspace(0,1/X.scalar(1),npoints);
else
    npoints=200*ceil(1/0.1);
    t1=linspace(0,1/0.1,npoints);
end
xt1=zeros(npoints,X.size_vector);

m=X.nodes;
for i=1:npoints
    t=t1(i);
    for j=1:X.size_vector
        xt1(i,j)=sum(X.vector(j,:).*exp(1i*[-m:m]*2*pi*t*X.scalar(1)));
    end
end
end


function plot_numerical_equilibria(lambda_hopf, equilibrium, dims)

% the lorenz map with appropriate parameters
f = @fn_Lorenz84; 
df = @derivatives_Lorenz84;

max_lambda = lambda_hopf + 0.007;
min_lambda = lambda_hopf - 0.007;

n_steps = 100;
steps = (max_lambda - min_lambda)/(n_steps-1);

vec_lambda = max_lambda:-steps:min_lambda;
lambda_iter = vec_lambda(1);
solutions = zeros(4,n_steps);
solutions(:,1) = Newton_handle(@(x)f(x,lambda_iter),equilibrium',@(x)df(x,lambda_iter));


for iter = 2:n_steps

    lambda_iter = vec_lambda(iter);
    solutions(:,iter) = Newton_handle(@(x)f(x,lambda_iter),...
        solutions(:,iter-1),@(x)df(x,lambda_iter));

end


plot3(vec_lambda,solutions(dims(1),:),solutions(dims(2),:),'c','LineWidth',1.5);
hold on
plot3(vec_lambda(n_steps/2),solutions(dims(1),n_steps/2),solutions(dims(2),n_steps/2),'r.','MarkerSize',30);


end
