function [xBar,iter, yBar,res,DFm,RES]=Newton_2(xBar,alpha,maxiter,min_res)
% function [xBar,yBar,res,DFm,RES]=Newton_Xi(xBar,alpha,maxiter,min_res)
% 
% INPUT:
% xBar      Xi_vector, first approximation of the solution in Fourier space;
% alpha     full_problem, that defines the system of ODEs considered;
% maxiter   positive integer, maximum number of iterations of the Newton
%           method (DEFAULT 20)
% min_res   positive real value, stopping criterion on the residual
%           (DEFAULT 10^-14)
%
% OUTPUT:
% xBar      Xi_vector, final solution;
% iter      number of iterations that were necessary
% yBar      Xi_vector, residual;
% res       positive real value, infinite norm of the residual.
% DFm       last derivative
% RES       sequence of residuals
%
% In the Newton method, the infinite norm on the Xi_vectors is used, that
% is, the biggest value between all is taken, considering scalars and
% elements of the Fourier transforms alike. 
% An error message is triggered if the derivative used is close to
% singularity. 

%global Display % if 1 display the final residual 
global use_intlab
global talkative

if talkative>2
    disp('Entering Newton')
end

if ~exist('maxiter','var') || isempty(maxiter)
    maxiter=10;
end
if ~exist('min_res','var') || isempty(min_res)
    min_res=10^-10;
end

if use_intlab
    warning('Newton method has poor convergence due to Intlab');
end

%disp('start')
RES=zeros(maxiter,1);
for iter=1:maxiter
    if talkative>3
        fprintf('Iteration %d, time %s\n',iter,datestr(now,13));
    end
    yBar=apply(alpha,xBar,0); % OK
    y=Xi_vec2vec(yBar); % working
    x=Xi_vec2vec(xBar);
    res=(norm(yBar));
    RES(iter,1:length(res))=res;
    res = max(res);
    if talkative>3
        fprintf('Residual %e, time %s\n',res,datestr(now,13));
    end
    %plot(i,res,'*');
    %hold on
    if res<min_res && iter>1
        break
    end
    
    DFm=derivative_to_matrix(derivative(alpha,xBar,0)); 
    if ~use_intlab && abs(det((DFm))) > 10e-5
        x=x-DFm\y;
    elseif use_intlab && abs(det(sup(DFm))) > 10e-5
        x=x-DFm\y;
    else
        error('derivative is singular, iteration %d',iter)
    end
    if any(isnan(x)) || any(isinf(x))
        error('Newton failed due to an infintiy or a NaN')
    end
    xBar=vec2Xi_vec(x,xBar.size_scalar,xBar.size_vector,xBar.nodes);
    xBar = symmetrise(xBar);
end

if iter==maxiter && res>min_res
    disp(res)
    figure
    semilogy(RES)
    error('Newton did not converge')
    return%
end

if talkative>2
    fprintf('The residual of the Newton method is %d.\n',res)
end
if talkative>2
    fprintf('Exiting Newton\n\n')
end

return