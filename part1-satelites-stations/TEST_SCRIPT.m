% Write to file:
fid = fopen('C:\Users\User\Desktop\hello.txt', 'a+');
a = ["a", "b", "c"];
fprintf(fid, '%i,%i,%i\n', 1:3);
fclose(fid);

% Read from file
% fid = fopen('C:\Users\User\Desktop\hello.txt', 'r');
% A = fscanf(fid,'%s');
% disp(A)
% fclose(fid);
% Export R table into latex table:
% print(xtable(newobject2, type = "latex"), file = "filename2.tex") 