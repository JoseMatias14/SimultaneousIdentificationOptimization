clear
close
clc

%Previously plant surface
load('PlantSurface');

% results
load('baseline');

% to control figure printing
% false - do not print
% true - print
pPDF = false;
pTIFF = false;
pr = [pPDF, pTIFF];

%% Plotting
inputs = {'u_{choke}','u_{comp}'};
measurements = {'Choke outlet pressure','Reservoir outflow','Compressor inflow','Compressor temperature','Compressor outlet pressure','Compressor efficiency','Compressor head','Compressor Power','Pump inflow','Pump power'};
units = {'P_{choke} [Pa]','q_{res} [m3/s]','q_{comp} [m3/s]','T_{comp} [K]','P_{comp} [Pa]','\nu [-]','H [m]','P_{ow} [W]','q_{pump} [m3/s]','P_{op} [W]'};
scaling = [1e7,1e0,1e0,1e2,1e7,1e0,1e1,1e4,1e-1,1e3];

%% 1. Input Comparison
leg = {'MAy','Opt'};

opt = load('PlantOptProfiles');

f1 = figure(1);

    subplot(2,1,1,'FontSize',10)
        plot(1:tEnd + 1,uOptArray(1,:),'k',1:tEnd,opt.uSP(1,1:tEnd),'k:')
        ylabel('u_{choke}','FontSize',10)
        xlabel('iteration','FontSize',10)
        xlim([1,tEnd + 1])
        %ylim([0.7,1])
    subplot(2,1,2,'FontSize',10)
        plot(1:tEnd + 1,uOptArray(2,:),'b',1:tEnd,opt.uSP(2,1:tEnd),'b:')
        ylabel('u_{comp}','FontSize',10)
        xlabel('iteration','FontSize',10) 
        xlim([1,tEnd + 1])
        %ylim([0.7,1])

        legend(leg,'Location','best','FontSize',9)

if pr(1)
    print(f1,'Results_Inputs','-dpdf')
end
if pr(2)
    print(f1,'-r1200','-dtiff','Results_Inputs.tif');
end


%% 2. comparing the actual gradient value with the estimated
leg = {'plant','hat'};

f2 = figure(2);
    for i = 1:ma.nMeas
        subplot(4,2,i,'FontSize',10)
            plot(1:tEnd,gradMeasPlantArray{i}(1,:),'k:')
            hold on
            plot(1:tEnd,gradPlantHatArray{i}(1,:),'k')

        ylabel(['g ',num2str(i)],'FontSize',10)

    end
    title('grad u choke')
    legend(leg,'Location','best','FontSize',9)

if pr(1)
    print(f2,'Results_Gradients_u1','-dpdf')
end
if pr(2)
    print(f2,'-r1200','-dtiff','Results_Gradients_u1.tif');
end

f3 = figure(3);
    for i = 1:ma.nMeas
        subplot(4,2,i,'FontSize',10)
            plot(1:tEnd,gradMeasPlantArray{i}(2,:),'b:')
            hold on
            plot(1:tEnd,gradPlantHatArray{i}(2,:),'b')

        ylabel(['g ',num2str(i)],'FontSize',10)

    end
title('grad u comp')
legend(leg,'Location','best','FontSize',9)

if pr(1)
    print(f3,'Results_Gradients_u2','-dpdf')
end
if pr(2)
    print(f3,'-r1200','-dtiff','Results_Gradients_u2.tif');
end

%% 3. comparing choosen model with actual model
leg = {'plant','hat'};
mod = {'H','D_1','D_2'};

%real plant array
rpArray = [];
for count = 1:tEnd
    if count < par.ph
        rpArray = [rpArray,1];
    elseif count < par.pd1 && count > par.ph_2_d1
        rpArray = [rpArray,2];
    elseif count > par.pd1_2_d2 
        rpArray = [rpArray,3];
    else 
        rpArray = [rpArray,0];
    end
end

f4 = figure(4);
subplot(2,1,1)
    stairs(1:tEnd,rpArray,'k:')
    hold on
    stairs(1:tEnd,modelArrayProb(2:end),'k')

    ylabel('Model','FontSize',10)
    title('Model Choice')
    legend(leg,'Location','best','FontSize',9)

 subplot(2,1,2)
    stairs(1:tEnd,maintenanceArray,'ko','markersize',3)

    xlim([0, 250])
    ylim([-0.5, 2.5])
    yticks([0,1,2])
    yticklabels({'Ok','Inspection','Maintenance'})
    ylabel('Model','FontSize',10)
    xlabel('SS Periods [-]','FontSize',10)
    title('Maintenance','FontSize',10)
hold off

if pr(1)
    print(f4,'Results_ModelChoice','-dpdf')
end
if pr(2)
    print(f4,'-r1200','-dtiff','Results_ModelChoice.tif');
end



% %% 4. comparing estimated parameters
% 
% %real plant array
% rpThetaArray{1} = [];
% rpThetaArray{2} = [];
% rpThetaArray{3} = [];
% 
% for count = 1:tEnd
%     if count < 31 
%         rpThetaArray{1} = [rpThetaArray{1},[0.582;-2.398;2.75;-3.969;4.303]];
%     elseif count < 91 && count > 60
%         rpThetaArray{2} = [rpThetaArray{2},[-0.5;1.8;-0.85]];
%     elseif count > 121
%         rpThetaArray{3} = [rpThetaArray{3},[0.43;0.7;5;1;2]];
%     else 
%         rpThetaArray{1} = [rpThetaArray{1},[0;0;0;0;0]];
%         rpThetaArray{2} = [rpThetaArray{2},[0;0;0]];
%         rpThetaArray{3} = [rpThetaArray{3},[0;0;0;0;0]];
%     end
% end
% 
% f5 = figure(5);
% color = ['k','b','r'];
% for i = 1:ma.nModels
%         plot(1:count,rpThetaArray{i},color(i),'LineSpec',':')
%         hold on
%         plot(1:count,thetaHat{i},color(i))
% 
% end
% 
% ylabel('Parameters','FontSize',10)
% title('Estimated parameters')
% print(f4,'-r80','-dpdf','ResultsParameters.pdf');
% 
% 
% 
