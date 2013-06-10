function cg_dumpfile(fn,txt)
% CG_DUMPFILE(FN,TXT) Write TXT to newly generated file FN.

fprintf(1,['Writing file: ',fn,'...']);
dlmwrite(fn,txt{1},'newline','pc','delimiter',''); % creates new file
for i = 2:length(txt)
    dlmwrite(fn,txt{i},'-append','newline','pc','delimiter','');
end
fprintf(1,'\tdone\n');