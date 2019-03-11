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
step_bound = 1000;
hold on; % keep plotting on the existing figure
title("Satelites orbiting the earth and some stations placed on the planet 's surface.");
earth_radius = 300;
axis([-earth_radius-200 earth_radius+200 -earth_radius-200 earth_radius+200]); %setting up figure size

% Creating satelite and station objects: remember, both take inputs as (arg_vel, arg_alt, arg_steps)
sat1 = satelite(30,earth_radius+60,8);
sat2 = satelite(-26,earth_radius+30,19);
stat1 = station(2,earth_radius,10);
stat2 = station(2,earth_radius,81);%for some reason, both stations overlap each other instead of being initialized at different positions!!!

for step=1:0.1:step_bound
    earth = viscircles([0 0],earth_radius,'Color','b'); %earth
    
    sat1.next_pos(step);
    sat2.next_pos(step);
    stat1.next_pos(step);
    stat2.next_pos(step);
   
    scoordinates1 = [sat1.get_long() sat1.get_lat()];
    scoordinates2 = [sat2.get_long() sat2.get_lat()];
    
    tcoordinates1 = [stat1.get_long() stat1.get_lat()];
    tcoordinates2 = [stat2.get_long() stat2.get_lat()];
   
    disp(tcoordinates1)%debug
    disp(tcoordinates2)%debug
    disp('--------------------------')%debug
   
    vis_sat1 = viscircles(scoordinates1,8,'Color','g');
    vis_sat2 = viscircles(scoordinates2,8,'Color','g');
    vis_stat1 = viscircles(tcoordinates1,12,'Color','r');
    vis_stat2 = viscircles(tcoordinates2,12,'Color','r');
       
    pause(0.02) %WARNING: all 'delete' (for graphic objects) functions must be after the 'pause' function 
    delete(vis_sat1) 
    delete(vis_sat2) 
    delete(vis_stat1) 
    delete(vis_stat2)
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
