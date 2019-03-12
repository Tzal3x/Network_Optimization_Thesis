%% Main file, here everything is assembled.
%{
Usefull hotkeys:
Ctrl+/ auto-comments the selected lines
Ctrl+T uncomments
%}

%{ 
https://uk.mathworks.com/help/optim/ug/fmincon.html#busog7r-2

fmincon will be out tool for solving this optimization problem
Syntax of fmincon: ------------------------------------------------------
x = fmincon(fun,x0,A,b)
x = fmincon(fun,x0,A,b,Aeq,beq)
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub)
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon)
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
x = fmincon(problem)
[x,fval] = fmincon(___)
[x,fval,exitflag,output] = fmincon(___)
[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(___)

[x,fval,exitflag,output,lambda,grad,hessian] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)

options = optimoptions('fmincon','Display','iter','Algorithm','sqp');
%}

%% Experimenting with fmincon
% fun = @(x)100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
% x0 = [-1,2];
% A = [1,2];
% b = 1;
% %constraint is in the form Ax=b
% x = fmincon(fun,x0,A,b) %#ok<*NOPTS>


%% Space set up:
hold on; % keep plotting on the existing figure
earth_radius = 300;
axis([-earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200]); %setting up figure size
title("Satelites orbiting the earth and some stations placed on the planet 's surface.");
step_bound = 50000;%should be the same as the 

sat1 = satelite(30,earth_radius+60,30,10,step_bound); %(arg_vel,arg_alt,arg_init_pos,arg_periods,arg_bound)
sat2 = satelite(60,earth_radius+80,65,10,step_bound);
sat3 = satelite(30,earth_radius+90,60,10,step_bound);
sat4 = satelite(80,earth_radius+50,55,10,step_bound);

station1 = satelite(10,earth_radius,90,10,step_bound);
station2 = satelite(10,earth_radius,10,10,step_bound);

for i=1:step_bound
    earth = viscircles([0 0],earth_radius,'Color','b'); %earth
    
    vis_sat1 = viscircles(sat1.lifetime_coordinates(:,i)',8,'Color','g');
    vis_sat2 = viscircles(sat2.lifetime_coordinates(:,i)',8,'Color','g');
    vis_sat3 = viscircles(sat3.lifetime_coordinates(:,i)',8,'Color','g');
    vis_sat4 = viscircles(sat4.lifetime_coordinates(:,i)',8,'Color','g');

    vis_station1 = viscircles(station1.lifetime_coordinates(:,i)',12,'Color','r');
    vis_station2 = viscircles(station2.lifetime_coordinates(:,i)',12,'Color','r');
    
    pause(0.05) %WARNING: all 'delete' (for graphic objects) functions must be after the 'pause' function 
    delete(vis_sat1) 
    delete(vis_sat2)
    delete(vis_sat3)
    delete(vis_sat4)
    
    delete(vis_station1)
    delete(vis_station2)
    
    delete(earth)
end
hold off;

%{
--------------------------------------------------------------------
Function definitions in a script must appear at the end of the file.
Only the first fuction can be accessed by other files. The rest are
acounted as local functions. MATLAB does no support default values on 
function parameters.
--------------------------------------------------------------------
%}
