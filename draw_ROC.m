function [] = draw_ROC(FPR_DIP_UP,TPR_DIP_UP,RECALL_DIP_UP,PRECISION_DIP_UP)
x=0:0.1:1;
figure(1)
plot(FPR_DIP_UP,TPR_DIP_UP,'r','LineWidth',1.6);

xlabel('FPR_DIP_UP');ylabel('TPR_DIP_UP');
title('ROC curves');
grid on;
axis on;
axis([0,1.0,0,1.0]);


hold on

x=0:0.1:1;
figure(2)
plot(RECALL_DIP_UP,PRECISION_DIP_UP,'g','LineWidth',1.6);
xlabel('RECALL_DIP_UP');ylabel('PRECISION_DIP_UP');%
title('PRÇúÏßÍ¼');
grid on;
axis on;
axis([0,1.0,0,0.5]);

% calculation of AUC and AUPR
AUC = trapz(FPR_DIP_UP,TPR_DIP_UP);
fprintf('\nAUC=%f\n\n',AUC);


AUPR = trapz(RECALL_DIP_UP,PRECISION_DIP_UP);
fprintf('AUPR=%f\n\n',AUPR);
end


