%This code run all numerical examples in this package

for cnt = 1:6
    fprintf('Begin example %d...\n',cnt);
    varSST(cnt,1);
    fprintf('Finish example %d, press enter to continue.\n',cnt);
    pause;
    close all;
end

fprintf('Begin noisy example ...\n');
varSSTns(1,1);
pause;
fprintf('Finish noisy example, press enter to continue.\n');
close all;

varSSTTwoStep(1,1);
