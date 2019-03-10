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


%% Test driving 
% n=600;
% x = 1:n;
% y = 1:n;      
% % x2 = 1:n;
% % y2 = 1:n;      
% z = repmat(408,length(x));
% sat = satelite(1,2,3);
% % sat2 = satelite(5,6,3);
% orbit = viscircles([sat.get_lat() sat.get_long()],408);
% for i = 1:n % 1*step, 2*step, 3*step (...)
%    sat.next_pos(i)
%    x(i) = sat.get_lat();
%    y(i) = sat.get_long();
%    disp([y(i),x(i)])
% %    sat2.next_pos(i)
% %    x2(i) = sat2.get_lat();
% %    y2(i) = sat2.get_long();
% end
% comet3(x,y,z) %comet (x,y) for 2d graph

%% Main start
step_bound = 1000;
hold on % keep plotting on the existing figure
title('Satelites orbiting the earth and some stations fixed in place')
rad = 408;
axis([-rad-100 rad+100 -rad-100 rad+100]) %setting up figure size
viscircles([0 0],408-88) %earth

% creating satelites:
sat1 = satelite(0,0,50);
sat2 = satelite(0,0,-60);
sat3 = satelite(0,0,-70);
sat4 = satelite(0,0,80);
sat5 = satelite(0,0,60);

sats = [sat1,sat2,sat3,sat4,sat5]; %list containing every satelite

for step=1:0.1:step_bound
    %Move every satelite:
    for sati = sats
        sati.next_pos(step);
        plot(sati.get_lat, sati.get_long,'r*');
    end
    drawnow
    pause(1) % 1 frame per second for animation purposes
end
hold off




%{
--------------------------------------------------------------------
Function definitions in a script must appear at the end of the file.
Only the first fuction can be accessed by other files. The rest are
acounted as local functions. MATLAB does no support default values on 
function parameters.
--------------------------------------------------------------------
%}