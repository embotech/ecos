function rowprint(x,name)

fprintf('%s =\n',name);
for i = 1:length(x)
    fprintf('%14.12e  ',x(i));
end
fprintf('\n');