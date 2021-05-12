clear all
close all
clc
load('results_standard')
% 

% nIter = 50;
% par = ParametersBlockReactor;
% ma.nModels = par.nr; %number of models

% %For plotting - previously calculated plant surface and optimum
load('PlantSurface_BR');


%% 1. ploting the results in plant surface
f1 = figure(1);
mod = {'M_{11}','M_{12}','M_{13}','M_{21}','M_{22}','M_{23}','M_{31}','M_{32}','M_{33}'};

iter = 1:nIter;

%plotting input profile
subplot(2,1,1)

cc = 0;
for ii = 1:size(modelArray,2)
    if modelArray(ii)~=4
        plot(cc,modelArray(ii),'k.','MarkerSize',10);
        hold on  
        %plot(cc,modelArray(ii),'kx','MarkerSize',8);
        cc = cc + 1;
    end
    
end

ax = gca;
ax.YGrid = 'on';
ax.GridLineStyle = '-';

yticks(1:ma.nModels)
yticklabels(mod)
ylim([0.8, ma.nModels+0.2])

xlabel('SS periods [-]','FontSize',10)
%ylabel('Model','FontSize',10)
title('Model with highest probability (M_{13} is correct)','FontSize',10)

%plotting input profile
subplot(2,1,2)
    %plot optimum

    area(0:nIter,modelArrayProbability')
    
    xlim([0, nIter])
    yticks(0:0.25:1)
    ylim([0 1])  
    
    xlabel('SS periods [-]')
    title('Probability of model')
    legend(mod,'Location','best','FontSize',6)

print(f1,'-r80','-dpdf','ResultsModels_stan.pdf');
hold off

%% 4. comparing the actual gradient value with the estimated

f4 = figure(4);
for i = 2
    plot(2:nIter + 1,gradPlantHatArray{i},'k')
    hold on
    plot(2:nIter + 1,gradMeasPlantArray{i},'k-.')
    
end

xticks(2:nIter + 1)
xlim([2 nIter + 1])

title('Objective Function Gradient','FontSize',10)
xlabel('SS periods [-]','FontSize',10)
ylabel('g Cc [bar/s]','FontSize',10)

leg = {'Plant','Estimation'};
legend(leg,'Location','best','FontSize',9)

% print(f7,'-r80','-dpdf','ResultsGradients.pdf');
hold off

%% 5. ploting the opt. tracking
f5 = figure(5);

subplot(2,1,2)
    plot(0:nIter,cost_OF*ones(1,nIter + 1),'k:','LineWidth',1.5);
    hold on
    plot(0:nIter,yPlantArray(2,1:3:end),'k','LineWidth',1.5);
    plot(0:nIter,yNoNoiseArray(2,1:3:end),'ok','MarkerSize',5);

    
    xlim([0, nIter-1])
    xticks(0:5:nIter-1)
    ylim([0.4 0.6])
    xlabel('SS periods [-]','FontSize',10)
    ylabel('C_C [mol/L]','FontSize',10)
    title('Objective function: \phi = C_C')
    
    leg = {'Optimal','Measured','Actual'};
    legend(leg,'Location','southeast','FontSize',8)
    
subplot(2,1,1)
    plot(0:nIter,uOptk*ones(1,nIter + 1),'k:','LineWidth',1.5);
    hold on
    stairs(0:nIter,inputPlantArray(1:3:end),'k','LineWidth',1.5);
    
    xlim([0, nIter-1])
    xticks(0:5:nIter-1)
    ylim([4 10])
    xlabel('SS periods [-]','FontSize',10)
    ylabel('F_1 [L/min]','FontSize',10)
    title('Manipulated variable: u = F_1')
    
    leg = {'Optimal','Standard'};
    legend(leg,'Location','northeast','FontSize',8)
    
print(f5,'-r80','-dpdf','ResultsOptTracking_stan.pdf');