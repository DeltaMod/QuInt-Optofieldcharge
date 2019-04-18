%% This function cleans up all current variables, and ensures that they are not greater than N in length%%
NClean = who;
Notthis = find(ismember(NClean, 'CONSOLEHISTORY'));
NClean(Notthis) = [];
fprintf(['Cleared: \n'])
if CleanWhat == "N is greater than Nold"
for i = 1:length(NClean)
    if length(eval(NClean{i})) > N
        fprintf(['       -', NClean{i},'\n'])
        clear (NClean{i})
    end
end
elseif CleanWhat == "N is smaller than Nold"
for i = 1:length(NClean)
    if length(eval(NClean{i})) < N && length(eval(NClean{i})) >= Nold
        fprintf(['       -', NClean{i},'\n'])
        clear (NClean{i})
    end
end
end
clear CleanWhat
clear NClean
clear Notthis
clear i