clear
close all 
clc

bl = load('baseline.mat');
g  = load('greedy.mat');
d  = load('dual.mat');
d005  = load('dual_005.mat');
d095  = load('dual_095.mat');

% to control figure printing
% false - do not print
% true - print
pPDF = true;
pTIFF = false;
pr = [pPDF, pTIFF];

%% Plotting
inputs = {'u_{choke}','u_{comp}'};
measurements = {'Choke outlet pressure','Reservoir outflow','Compressor inflow','Compressor temperature','Compressor outlet pressure','Compressor efficiency','Compressor head','Compressor Power','Pump inflow','Pump power'};
units = {'P_{choke} [Pa]','q_{res} [m3/s]','q_{comp} [m3/s]','T_{comp} [K]','P_{comp} [Pa]','\nu [-]','H [m]','P_{ow} [W]','q_{pump} [m3/s]','P_{op} [W]'};
scaling = [1e7,1e0,1e0,1e2,1e7,1e0,1e1,1e4,1e-1,1e3];

%% 1. Input Comparison (separate)
leg = {'Baseline','Greedy','Dual'};

f1 = figure(1);

    subplot(3,1,1,'FontSize',10)
        plot(1:bl.tEnd + 1,bl.uOptArray(1,:),'bo','MarkerSize',2);
        title('Baseline')
        ylabel('u_{choke} [%]','FontSize',10)
        xlabel('SS periods [-]','FontSize',10)
        xlim([1,bl.tEnd + 1])
        %ylim([0.7,1])
        
    subplot(3,1,2,'FontSize',10)
        plot(1:g.tEnd + 1,g.uOptArray(1,:),'-r*','MarkerSize',2);
        title('MAP greedy')
        ylabel('u_{choke} [%]','FontSize',10)
        xlabel('SS periods [-]','FontSize',10)
        xlim([1,g.tEnd + 1])
        %ylim([0.7,1])
        
    subplot(3,1,3,'FontSize',10)
        plot(1:d.tEnd + 1,d.uOptArray(1,:),'-kd','MarkerSize',2);
        title('Dual Control')
        ylabel('u_{comp} [%]','FontSize',10)
        xlabel('SS periods [-]','FontSize',10)
        xlim([1,d.tEnd + 1])
        %ylim([0.7,1])

if pr(1)
    print(f1,'Results_Inputs','-dpdf')
end
if pr(2)
    print(f1,'-r1200','-dtiff','Results_Inputs.tif');
end

%% 1.1 Input Comparison (together)
f2 = figure(2);

        plot(1:bl.tEnd + 1,bl.uOptArray(1,:),'-bo','MarkerSize',2);
        hold on
        plot(1:g.tEnd + 1,g.uOptArray(1,:),'-r*','MarkerSize',2);
        plot(1:d.tEnd + 1,d.uOptArray(1,:),'-kd','MarkerSize',2);
        
        title('Input Comparison')
        ylabel('u_{choke} [%]','FontSize',10)
        xlabel('SS periods [-]','FontSize',10)
        xlim([1,bl.tEnd + 1])
        legend(leg,'Location','best','FontSize',9)

if pr(1)
    print(f2,'Results_Inputs2','-dpdf')
end
if pr(2)
    print(f2,'-r1200','-dtiff','Results_Inputs2.tif');
end

%% 1.2 Input Comparison (Dual)
f3 = figure(3);
leg2 = {'Dual','\alfa = 0.05','\alfa = 0.95'};

        plot(1:d.tEnd + 1,d.uOptArray(1,:),'kd','MarkerSize',3);
        hold on
        plot(1:d005.tEnd + 1,d005.uOptArray(1,:),':b');
       	plot(1:d095.tEnd + 1,d095.uOptArray(1,:),'--b','MarkerSize',2);

        title('Comparison \alpha')
        ylabel('u_{choke} [%]','FontSize',10)
        xlabel('SS periods [-]','FontSize',10)
        xlim([1,d.tEnd + 1])
        legend(leg2,'Location','best','FontSize',9)

if pr(1)
    print(f3,'Results_Inputs3','-dpdf')
end
if pr(2)
    print(f3,'-r1200','-dtiff','Results_Inputs3.tif');
end

%% 2. Probability vs. Chosen Model
f4 = figure(4);
leg3 = {'H','D_1','D_2'};

%real plant array
rpArray = [];
for count = 1:bl.tEnd
    if count < bl.par.ph
        rpArray = [rpArray,1];
    elseif count < bl.par.pd1 && count > bl.par.ph_2_d1
        rpArray = [rpArray,2];
    elseif count > bl.par.pd1_2_d2 
        rpArray = [rpArray,3];
    else 
        rpArray = [rpArray,0];
    end
end

subplot(3,1,1)
    stairs(1:bl.tEnd,rpArray,'r:','LineWidth',1.5)
    hold on
    stairs(1:bl.tEnd,bl.modelArrayProb(2:end),'kx','MarkerSize',4)

    xlim([0, g.tEnd])
    ylim([0.5, 3.5])
    yticks([1,2,3])
    yticklabels({'H','D1','D2'})
    ylabel('Model','FontSize',10)
    title('Model Choice (Baseline)')

%plotting probability profile
subplot(3,1,2)
    area(0:g.tEnd,g.probModelArray')
    
    xlim([0, g.tEnd])
    ylim([0 1])  
    
    title('Probability of model (Greedy)')
    xlabel('SS periods [-]')
    legend(leg3,'FontSize',6,'Location','best')
    
    
%plotting probability profile
subplot(3,1,3)
    area(0:d.tEnd,d.probModelArray')
    
    xlim([0, d.tEnd])
    ylim([0 1])  
    
    title('Model probability (Dual)')
    xlabel('SS periods [-]')
    ylabel('Probability of model')
    legend(leg3,'FontSize',6,'Location','best')

if pr(1)
    print(f4,'Results_ModelChoice','-dpdf')
end
if pr(2)
    print(f4,'-r1200','-dtiff','Results_ModelChoice.tif');
end

%% 3. Profit Comparison (together)

f5 = figure(5);
    %computing the value
    OF_b = bl.yPlantArray(3,1:3:end)./(bl.yPlantArray(6,1:3:end));
    OF_g = g.yPlantArray(3,1:3:end)./(g.yPlantArray(6,1:3:end));
    OF_d = d.yPlantArray(3,1:3:end)./(d.yPlantArray(6,1:3:end));

    plot(1:bl.tEnd + 1,OF_b,'-bo','MarkerSize',2);
    hold on
    plot(1:g.tEnd + 1,OF_g,'-r*','MarkerSize',2);
    plot(1:d.tEnd + 1,OF_d,'-kd','MarkerSize',2);

    title('Comparison (Economic OF)')
    ylabel('J','FontSize',10)
    xlabel('SS periods [-]','FontSize',10)
    xlim([1,bl.tEnd + 1])
    legend(leg,'Location','best','FontSize',9)

if pr(1)
    print(f5,'Results_Profit','-dpdf')
end
if pr(2)
    print(f5,'-r1200','-dtiff','Results_Profit.tif');
end

%% 3. Profit Comparison (alone)

f6 = figure(6);

subplot(3,1,1), 
    plot(1:g.tEnd + 1,OF_g - OF_b,'-kx','MarkerSize',2)
    grid on
    title('OF dif. (Greedy - Baseline)')
    ylabel('J','FontSize',10)
    xlabel('SS periods [-]','FontSize',10)
    xlim([1,g.tEnd + 1])
    
subplot(3,1,2), 
    plot(1:d.tEnd + 1,OF_d - OF_b,'-kx','MarkerSize',2)
    grid on
    title('OF dif. (Dual - Baseline)')
    ylabel('J','FontSize',10)
    xlabel('SS periods [-]','FontSize',10)
    xlim([1,d.tEnd + 1])
    
subplot(3,1,3)
    plot(1:g.tEnd + 1,OF_g - OF_d,'-kx','MarkerSize',2)
    grid on
    title('OF dif.(Greedy - Dual)')
    ylabel('J','FontSize',10)
    xlabel('SS periods [-]','FontSize',10)
    xlim([1,g.tEnd + 1])
    
if pr(1)
    print(f6,'Results_Profit','-dpdf')
end
if pr(2)
    print(f6,'-r1200','-dtiff','Results_Profit.tif');
end