sep = [1,0.5,0.1,0.05,0.01,0.005];
load('gammaConditionTest.mat')

fill([-log(sep),fliplr(-log(sep))],[log(condK_gammaP1), fliplr(log(condK_gamma10))],'r','FaceAlpha',0.3)
hold on
fill([-log(sep),fliplr(-log(sep))],[log(condLKL_gammaP1), fliplr(log(condLKL_gamma10))],'g','FaceAlpha',0.3)
ylabel('log(condition #)','fontsize',20)
xlabel('-log(sep)','fontsize',20)