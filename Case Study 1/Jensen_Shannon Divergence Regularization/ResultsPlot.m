%% Uncomment if using the script directly (without running Main.m)
% clear all
% close all
% clc
load('results_Dual') 
% par = ParametersBlockReactor;
% ma.nModels = par.nr; %number of models
% load('PlantSurface_BR');

% %For plotting - previously calculated plant surface and optimum
load('PlantSurface_BR');

%%
% to control figure printing
% false - do not print
% true - print
pPDF = true;
pTIFF = false;
pr = [pPDF, pTIFF];

%% 1. plotting model choice and associated probabilities
f1 = figure(1);
mod = {'M_{11}','M_{12}','M_{13}','M_{21}','M_{22}','M_{23}','M_{31}','M_{32}','M_{33}'};

iter = 1:nIter;

%plotting input profile
subplot(2,1,1)

cc = 0;
for ii = 1:size(modelArrayProb,2) 
    if modelArrayProb(ii)~=10 %probing points     
        plot(cc,modelArrayProb(ii),'k.','MarkerSize',10);
        hold on
        
        cc = cc + 1;
    end
    
end

cc = 0;
% for ii = 1:size(probModelArray,2) 
%     [~,maxTemp] = max(probModelArray(:,ii));
%     plot(cc,maxTemp,'kx','MarkerSize',10);
%     hold on
%     
%     cc = cc + 1;
% end

ax = gca;
ax.YGrid = 'on';
ax.GridLineStyle = '-';

yticks(1:ma.nModels)
yticklabels(mod)
ylim([0.8, ma.nModels+0.2])
%ylabel('Model','FontSize',10)

xlabel('SS periods [-]','FontSize',10)
title('Model with highest probability (M_{13} is correct)','FontSize',10)

%plotting probability profile
subplot(2,1,2)
    area(0:nIter,probModelArray')
    
    xlim([0, nIter])
    yticks(0:0.25:1)
    ylim([0 1])  
    
    xlabel('SS periods [-]')
    title('Probability of model')
    legend(mod,'FontSize',6,'Location','best')

hold off

if pr(1)
    print(f1,'-r80','-dpdf','ResultsModels_DUAL.pdf');
elseif pr(2)
    print(f1,'-r1200','-dtiff','ResultsModels_DUAL.tif');
end

%% 2. comparing the actual gradient value with the estimated

f4 = figure(4);
for i = 2
    plot(2:nIter+1,gradPlantHatArray{i},'k')
    hold on
    plot(2:nIter+1,gradMeasPlantArray{i},'k-.')
    
end

xticks(2:nIter+1)
xlim([2 nIter+1])
ylim([-0.5 0.5]) 

title('Objective Function Gradient','FontSize',10)
xlabel('SS periods [-]','FontSize',10)
ylabel('g Cc [bar/s]','FontSize',10)

leg = {'Plant','Estimation'};
legend(leg,'Location','best','FontSize',9)

hold off
% if pr(1)
%     print(f4,'-r80','-dpdf','ResultsGradients.pdf');
% elseif pr(2)
%     print(f4,'-r1200','-dtiff','ResultsGradients.tif');
% end


%% 3. plotting opt. tracking
f5 = figure(5);

subplot(2,1,2)
    plot(0:nIter,cost_OF*ones(1,nIter+1),'k:','LineWidth',1.5);
    hold on
    plot(0:nIter,yPlantArray(2,1:3:end),'r','LineWidth',1.5);
    plot(0:nIter,yNoNoiseArray(2,1:3:end),'or','MarkerSize',5);

    xlim([0, nIter])
    xticks(0:5:nIter)
    ylim([0.4 0.6])
    xlabel('SS periods [-]','FontSize',10)
    ylabel('C_C [mol/L]','FontSize',10)
    title('Objective function: \phi = C_C')
    
    leg = {'Optimal','Measured','Actual'};
    legend(leg,'Location','southeast','FontSize',8)
    
subplot(2,1,1)
    h = area((1:nIter)',[min(ukLambdaArray,[],1)',max(ukLambdaArray,[],1)' - min(ukLambdaArray,[],1)'],'LineStyle','none');
    h(1).FaceColor = [1 1 1];
    h(1).HandleVisibility = 'off';
    h(2).HandleVisibility = 'off';
    h(2).FaceColor = [1 0 0];
    h(2).FaceAlpha = 0.2;
    hold on

    plot(0:nIter,uOptk*ones(1,nIter+1),'k:','LineWidth',1.5);
    hold on
    stairs(0:nIter,inputPlantArray(1:3:end),'rx','LineWidth',1.5);
    st = load('results_standard');
    stairs(0:nIter,st.inputPlantArray(1:3:end),'k','LineWidth',1.5);
        
    xlim([0, nIter])
    xticks(0:5:nIter)
    ylim([4 10])
    xlabel('SS periods [-]','FontSize',10)
    ylabel('F_1 [L/min]','FontSize',10)
    title('Manipulated variable: u = F_1')
    
    leg = {'Optimal','JSDR','Baseline'};  
    legend(leg,'Location','northeast','FontSize',9)

    
if pr(1)
    print(f5,'-r80','-dpdf','ResultsOptTracking_DUAL.pdf');
elseif pr(2)
    print(f5,'-r1200','-dtiff','ResultsOptTracking_DUAL.tif');
end   
    
%% Plotting the input regions     
figure(6)
    h = area((1:nIter)',[min(ukLambdaArray,[],1)',max(ukLambdaArray,[],1)' - min(ukLambdaArray,[],1)'],'LineStyle','none');
    h(1).FaceColor = [1 1 1];
    h(1).HandleVisibility = 'off';
    h(2).FaceColor = [0 0 0];
    h(2).FaceAlpha = 0.2;
    hold on
    
    scatter(1:nIter,inputPlantArray(4:3:end),40,'x','MarkerEdgeColor',[1, 0, 0],'LineWidth',1)
    scatter(1:nIter,ukLambdaArray(12,:),60,'o','MarkerEdgeColor',[0, 0.5, 0],'LineWidth',1)    
    scatter(1:nIter,ukLambdaArray(1,:),60,'o','MarkerEdgeColor',[0.75, 0.75, 0],'LineWidth',1) 
                                   
    plot(1:nIter,uOptk*ones(1,nIter),'k:','LineWidth',1.5);
    
    xlim([0, nIter])
    xticks(0:5:nIter)
    ylim([4.2 8.5])
    xlabel('SS periods [-]','FontSize',10)
    ylabel('F_1 [L/min]','FontSize',10)
    title('Manipulated variable: u = F_1')
    
    leg = {'Input Region','Implemented','J','D','Optimal'};
    legend(leg,'Location','northeast','FontSize',9)
 
if pr(1)
    print(f5,'-r80','-dpdf','ResultsOptRegion_DUAL.pdf');
elseif pr(2)
    print(f5,'-r1200','-dtiff','ResultsOptRegion_DUAL.tif');
end   
    