load('betaConditionTest.mat')

fill([-log(sep),fliplr(-log(sep))],[log(condK_beta5), fliplr(log(condK_betaP25))],'r','FaceAlpha',0.3)
hold on
fill([-log(sep),fliplr(-log(sep))],[log(condLKL_beta5), fliplr(log(condLKL_betaP25))],'g','FaceAlpha',0.3)
ylabel('log(condition #)','fontsize',20)
xlabel('-log(sep)','fontsize',20)