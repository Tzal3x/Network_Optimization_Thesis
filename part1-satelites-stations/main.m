%% Main file, here everything is ammased.
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
x = [];
y = [];
sat = satelite(1,2,3);

for i = 1:600
   sat.next_pos(i)
   x = [x,sat.get_lat()];
   y = [y,sat.get_long()];
%    pause(0.2) 
end

comet(x,y)


%{
--------------------------------------------------------------------
Function definitions in a script must appear at the end of the file.
Only the first fuction can be accessed by other filies. The rest are
acounted as local functions.
--------------------------------------------------------------------
%}
function out = f(x) %note: MATLAB does no support default values on function parameters
out = 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;
end