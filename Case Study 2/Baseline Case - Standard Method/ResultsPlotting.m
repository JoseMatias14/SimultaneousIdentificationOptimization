% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % For plotting without running the simulation, uncomment this section %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% clear
% close
% clc
% 
% % results
% load('results_standard');

%%
% loading precomputed optimal profiles
opt = load('PlantOptProfiles');

% to control figure printing
% false - do not print
% true - print
pPDF = true;
pTIFF = false;
pr = [pPDF, pTIFF];

%% names and tags
inputs = {'u_{choke}','u_{comp}'};
measurements = {'Choke outlet pressure','Reservoir outflow','Compressor inflow','Compressor temperature','Compressor outlet pressure','Compressor efficiency','Compressor head','Compressor Power','Pump inflow','Pump power'};
units = {'P_{choke} [Pa]','q_{res} [m3/s]','q_{comp} [m3/s]','T_{comp} [K]','P_{comp} [Pa]','\nu [-]','H [m]','P_{ow} [W]','q_{pump} [m3/s]','P_{op} [W]'};
scaling = [1e7,1e0,1e0,1e2,1e7,1e0,1e1,1e4,1e-1,1e3];
mod = {'H','D_1','D_2'};


%% 1. Input Comparison
leg = {'MAy','Opt'};

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
    print(f1,'Results_Inputs_STD','-dpdf')
end
if pr(2)
    print(f1,'-r1200','-dtiff','Results_Inputs_STD.tif');
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

% if pr(1)
%     print(f2,'Results_Gradients_u1','-dpdf')
% end
% if pr(2)
%     print(f2,'-r1200','-dtiff','Results_Gradients_u1.tif');
% end

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

% if pr(1)
%     print(f3,'Results_Gradients_u2','-dpdf')
% end
% if pr(2)
%     print(f3,'-r1200','-dtiff','Results_Gradients_u2.tif');
% end

%% 3. comparing choosen model with actual model
leg = {'plant','hat'};

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

hold off

if pr(1)
    print(f4,'Results_ModelChoice_STD','-dpdf')
end
if pr(2)
    print(f4,'-r1200','-dtiff','Results_ModelChoice_STD.tif');
end
