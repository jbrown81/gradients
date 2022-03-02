function [y_pos_all,y_vel_all]=gradient_ode(grad_slopes,grad_slope_deltas,all_betas,pred_interval)
% gradient_ode.m
% Jesse Brown
% 01/2022
% jesse.brown@ucsf.edu

%   [y_pos_all,y_vel_all]=gradient_ode(grad_slopes,grad_slope_deltas,all_betas,pred_interval)
%   returns the solution to the system of second order linear ordinary
%   differential equations, ie the forecasted 'position' (y_pos_all) and
%   'velocity' (y_vel_all). The solution is based on the initial conditions 
%   for gradient position (ie the gradient timeseries, grad_slopes) and velocity
%   (ie first derivatives, grad_slope_deltas), the coupling parameters between 
%   gradients (all_betas), and the time interval over which to make predictions (pred_interval).

stepsize=1;
tspan=[0:stepsize:pred_interval];
n_tps=size(grad_slopes,1);
n_comps=size(grad_slopes,2);

% y_pos_all is the ode solution
% aka the forecasted gradient timeseries from the initial conditions
y_pos_all=zeros(length(tspan),n_comps,n_tps);
count=1;
for i=1:n_tps
    t0=[];
    for j=1:n_comps
        t0=[t0;grad_slopes(i,j);grad_slope_deltas(i,j);];
    end
    
    % solve system of second order differential equations
    options=odeset('reltol',1e-5);
    [t,y]=ode45(@(t,y) f(t,y,all_betas),tspan,t0,options);
    
    y_pos=y(:,1:2:(n_comps*2)); % gradient positions
    y_veloc=y(:,2:2:(n_comps*2)); % gradient velocities
    y_pos_all(:,:,count)=y_pos;
    y_vel_all(:,:,count)=y_veloc;
    
    count=count+1;
    if ~mod(count,500)
        disp(count)
    end
end

function vals=f(t,y,all_betas)
% represent second order differential equation as system of two first order
% equations in vector format
% see https://www.mathworks.com/videos/solving-odes-in-matlab-8-systems-of-equations-117652.html
% do this for all components
% return values array that ode45 uses to compute numerical solution
n_comps=size(all_betas,1);
vals=[];
for i=2:2:n_comps*2
    cur_comp_num=i/2;
    cur_betas_int=all_betas(cur_comp_num,1);
    cur_betas_params=all_betas(cur_comp_num,2:end);
    % first value in 'vals' is x'', aka y(2)'
    vals=[vals;y(i);];
    cur_eq_val=cur_betas_int;
    for j=1:length(cur_betas_params)
        cur_eq_val=cur_eq_val+cur_betas_params(j)*y(j);
    end
    % second value in 'vals' is what x'' is equal to, aka the differential equation
    vals=[vals;cur_eq_val;];
end
end

end

