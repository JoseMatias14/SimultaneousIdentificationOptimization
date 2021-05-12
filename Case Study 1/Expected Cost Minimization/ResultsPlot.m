%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For plotting without running the simulation, uncomment this section %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all
% close all
% clc
% 
% % loading simulation results
% load('results_ECM')
% 
% nIter = 50;
% par = ParametersBlockReactor;
% ma.nModels = par.nr; %number of models

%%
%For plotting - previously calculated plant surface and optimum
load('PlantSurface_BR'); 

% Loading results from the baseline case solution for comparison
st = load('results_standard');

% to control figure printing
% false - do not print
% true - print
pPDF = true;
pTIFF = false;
pr = [pPDF, pTIFF];

%% names and tags
mod = {'M_{11}','M_{12}','M_{13}','M_{21}','M_{22}','M_{23}','M_{31}','M_{32}','M_{33}'};

%% 1. plotting model with the highest probability + evoution of posteriori belief
f1 = figure(1);

% vector with the SS periods
iter = 1:nIter;

% plotting model with the highest probability
subplot(2,1,1)

    % looping through the model probability vector and excluding the
    % probing points
    cc = 0;
    for ii = 1:size(modelArrayProb,2) 
        if modelArrayProb(ii)~=0  % excluding probing points   
            plot(cc,modelArrayProb(ii),'k.','MarkerSize',10);
            hold on

            cc = cc + 1;
        end

    end

    ax = gca;
    ax.YGrid = 'on';
    ax.GridLineStyle = '-';

    yticks(1:ma.nModels)
    yticklabels(mod)
    ylim([0.8, ma.nModels+0.2])
    %ylabel('Model','FontSize',10)

    xlabel('SS periods [-]','FontSize',10)
    title('Model with highest probability (M_{13} is correct)','FontSize',10)

%plotting probability profiles
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
    print(f1,'-r80','-dpdf','ResultsModels_ECM.pdf');
elseif pr(2)
    print(f1,'-r1200','-dtiff','ResultsModels_ECM.tif');
end

%% 2. comparing the actual gradient value with the estimated
f2 = figure(2);
for i = 2
    plot(2:nIter+1,gradPlantHatArray{i},'k')
    hold on
    plot(2:nIter+1,gradMeasPlantArray{i},'k-.')
    
end

xticks(2:4:nIter+1)
xlim([2 nIter+1])

title('Objective Function Gradient','FontSize',10)
xlabel('SS periods [-]','FontSize',10)
ylabel('g Cc [bar/s]','FontSize',10)

leg = {'Plant','Estimation'};
legend(leg,'Location','best','FontSize',9)

% if pr(1)
%     print(f2,'-r80','-dpdf','ResultsGradients.pdf');
% elseif pr(2)
%     print(f2,'-r1200','-dtiff','ResultsGradients.tif');
% end
hold off

%% 3. ploting the opt. tracking
f3 = figure(3);

% INPUTS
subplot(2,1,1)
    plot(0:nIter,uOptk*ones(1,nIter+1),'k:','LineWidth',1.5);
    hold on
    stairs(0:nIter,inputPlantArray(1:3:end),'b','LineWidth',1.5);
    
%     st = load('results_standard');
    stairs(0:nIter,st.inputPlantArray(1:3:end),'k','LineWidth',1.5);
    
    xlim([0, nIter])
    xticks(0:5:nIter)
    ylim([4 10])
    xlabel('SS periods [-]','FontSize',10)
    ylabel('F_1 [L/min]','FontSize',10)
    title('Manipulated variable: u = F_1')
    
    leg = {'Optimal','ECM','Baseline'};
    legend(leg,'Location','northeast','FontSize',8)

% OF: C concentration 
subplot(2,1,2)
    plot(0:nIter,cost_OF*ones(1,nIter+1),'k:','LineWidth',1.5);
    hold on
    plot(0:nIter,yPlantArray(2,1:3:end),'b','LineWidth',1.5);
    plot(0:nIter,yNoNoiseArray(2,1:3:end),'ob','MarkerSize',5);

    xlim([0, nIter])
    xticks(0:5:nIter)
    ylim([0.4 0.6])
    xlabel('SS periods [-]','FontSize',10)
    ylabel('C_C [mol/L]','FontSize',10)
    title('Objective function: \phi = C_C')
    
    leg = {'Optimal','Measured','Actual'};
    legend(leg,'Location','southeast','FontSize',8)
    
if pr(1)
    print(f3,'-r80','-dpdf','ResultsOptTracking_ECM.pdf');
elseif pr(2)
    print(f3,'-r1200','-dtiff','ResultsOptTracking_ECM.tif');
end