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

%% gambiarra para cortar os valores
tEnd = 220;

blThres = 161;
gThres = 161;
dThres = 212;

% inputs (all have the same size)
[temp1,temp2] = size(bl.uOptArray);
bl.uOptArray(:,blThres:end) = zeros(temp1,length(blThres:temp2));
g.uOptArray(:,gThres:end) = zeros(temp1,length(gThres:temp2));
d.uOptArray(:,dThres:end) = zeros(temp1,length(dThres:temp2));
d005.uOptArray(:,dThres:end) = zeros(temp1,length(dThres:temp2));
d095.uOptArray(:,dThres:end) = zeros(temp1,length(dThres:temp2));

% model probability (all have the same size)
[temp3,temp4] = size(bl.modelArrayProb);
bl.modelArrayProb(:,blThres:end) = zeros(temp3,length(blThres:temp4));
g.modelArrayProb(:,gThres:end) = zeros(temp3,length(gThres:temp4));
d.modelArrayProb(:,dThres:end) = zeros(temp3,length(dThres:temp4));
d005.modelArrayProb(:,dThres:end) = zeros(temp3,length(dThres:temp4));
d095.modelArrayProb(:,dThres:end) = zeros(temp3,length(dThres:temp4));

% measurements (all have the same size)
[temp5,temp6] = size(bl.yPlantArray);
bl.yPlantArray(:,3*blThres:end) = zeros(temp5,length(3*blThres:temp6));
g.yPlantArray(:,3*gThres:end) = zeros(temp5,length(3*gThres:temp6));
d.yPlantArray(:,3*dThres:end) = zeros(temp5,length(3*dThres:temp6));
d005.yPlantArray(:,3*dThres:end) = zeros(temp5,length(3*dThres:temp6));
d095.yPlantArray(:,3*dThres:end) = zeros(temp5,length(3*dThres:temp6));

%% 1. Input Comparison (separate)
leg = {'Baseline','Greedy','Dual'};

f1 = figure(1);

    subplot(3,1,1,'FontSize',10)
%         area([0, 50],[1, 1],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2,'LineStyle','none')
%         hold on
%         area([100, 150],[1, 1],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2,'LineStyle','none')
%         area([200, 250],[1, 1],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2,'LineStyle','none')
%         
       
        a = area([gThres, 220],[0.8 0.8],'FaceColor','w','LineWidth',1.5);
        a.FaceAlpha = 0.2;
        hold on 
        text(190,0.65,'Shutdown','HorizontalAlignment','center');

    
        stairs(1:bl.tEnd + 1,bl.uOptArray(1,:),'k','LineWidth',2.5);
        %
        hold on
        yline(0.7,'Color',[0.7,0.7,0.7]);
        yline(0.6,'Color',[0.7,0.7,0.7]);
        
        ax = gca;
        ax.YGrid = 'on';
        ax.GridLineStyle = '-';
        
        title('Manipulated variable (Baseline)')
        ylabel('v_{ch} [%]','FontSize',10)
        xlabel('SS periods [-]','FontSize',10)
        xlim([1,tEnd])
        ylim([0.5,0.8])
        yticks(0.5:0.1:0.8)
        
    subplot(3,1,2,'FontSize',10)
%         area([0, 50],[1, 1],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2,'LineStyle','none')
%         hold on
%         area([100, 150],[1, 1],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2,'LineStyle','none')
%         area([200, 250],[1, 1],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2,'LineStyle','none')

        a = area([gThres, 220],[0.8 0.8],'FaceColor','w','LineWidth',1.5);
        a.FaceAlpha = 0.2;
        hold on 
        text(190,0.65,'Shutdown','HorizontalAlignment','center');
        
        stairs(1:g.tEnd + 1,g.uOptArray(1,:),'b','LineWidth',2.5);
        %         
        hold on 
        yline(0.7,'Color',[0.7,0.7,0.7]);
        yline(0.6,'Color',[0.7,0.7,0.7]);
        
        ax = gca;
        ax.YGrid = 'on';
        ax.GridLineStyle = '-';
        
        title('Manipulated variable (ECM)')
        ylabel('v_{ch} [%]','FontSize',10)
        xlabel('SS periods [-]','FontSize',10)
        xlim([1,tEnd])
        ylim([0.5,0.8])
        yticks(0.5:0.1:0.8)
        
    subplot(3,1,3,'FontSize',10)
%         area([0, 50],[1, 1],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2,'LineStyle','none')
%         hold on
%         area([100, 150],[1, 1],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2,'LineStyle','none')
%         area([200, 250],[1, 1],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2,'LineStyle','none')

        a = area([dThres, 220],[0.8 0.8],'FaceColor','w','LineWidth',1.5);
        a.FaceAlpha = 0.2;
        hold on 
        text(215,0.65,'Shutdown','HorizontalAlignment','center','Rotation',90);

        h = area((1:d005.tEnd + 1)',[d095.uOptArray(1,:)',d005.uOptArray(1,:)' - d095.uOptArray(1,:)'],'LineStyle','none');
        h(1).FaceColor = [1 1 1];
        h(1).HandleVisibility = 'off';
        h(2).HandleVisibility = 'off';
        h(2).FaceColor = [1 0 0];
        h(2).FaceAlpha = 0.2;        
        hold on 
        yline(0.7,'Color',[0.7,0.7,0.7]);
        yline(0.6,'Color',[0.7,0.7,0.7]);
        
        
        stairs(1:d.tEnd + 1,d.uOptArray(1,:),'rx','LineWidth',0.5);
        
        title('Manipulated variable (JSDR)')
        ylabel('v_{ch} [%]','FontSize',10)
        xlabel('SS periods [-]','FontSize',10)
        xlim([1,tEnd])
        ylim([0.5,0.8])
        yticks(0.5:0.1:0.8)

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
        ylabel('v_{ch} [%]','FontSize',10)
        xlabel('SS periods [-]','FontSize',10)
        xlim([1,tEnd])
        ylim([0.45,0.85])
        legend(leg,'Location','best','FontSize',9)

if pr(1)
    print(f2,'Results_Inputs2','-dpdf')
end
if pr(2)
    print(f2,'-r1200','-dtiff','Results_Inputs2.tif');
end

%% 1.2 Input Comparison (Dual)
f3 = figure(3);
leg2 = {'Dual','\alpha = 0.05','\alpha = 0.95'};

        plot(1:d.tEnd + 1,d.uOptArray(1,:),'kd','MarkerSize',3);
        hold on
        plot(1:d005.tEnd + 1,d005.uOptArray(1,:),':b');
       	plot(1:d095.tEnd + 1,d095.uOptArray(1,:),'--b','MarkerSize',2);

        title('Comparison \alpha')
        ylabel('v_{ch} [%]','FontSize',10)
        xlabel('SS periods [-]','FontSize',10)
        xlim([1,tEnd])
        ylim([0.45,1])
        legend(leg2,'Location','best','FontSize',9)

if pr(1)
    print(f3,'Results_Inputs3','-dpdf')
end
if pr(2)
    print(f3,'-r1200','-dtiff','Results_Inputs3.tif');
end

%% 2. Probability vs. Chosen Model A
f4 = figure(4);
leg3 = {'H','D_1','D_2'};

subplot(3,1,1)
    %%%%%%%%%%
    % GAMBIS %
    %%%%%%%%%%
    %hold on
    a = area([gThres - 1, 220],[4.5 4.5],'FaceColor','w','LineWidth',1.5);
    a.FaceAlpha = 0.2;
    hold on 
    area([0, 50],[5, 5],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2,'LineStyle','none')
    area([100, 150],[5, 5],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2,'LineStyle','none')
    area([200, 250],[5, 5],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.1,'LineStyle','none')
    
    stairs(1:bl.tEnd,bl.modelArrayProb(2:end),'k.','MarkerSize',10)
    
    text(25,4,'H','HorizontalAlignment','center')
    text(125,4,'D_1','HorizontalAlignment','center')
    text(210,4,'D_2','HorizontalAlignment','center','Color',[0.5,0.5,0.5])
    text(190,2.5,'Shutdown','HorizontalAlignment','center')
    
    ax = gca;
    ax.YGrid = 'on';
    ax.GridLineStyle = '-';
    
    xlim([0, tEnd])
    ylim([0.5, 4.5])
    yticks([1,2,3])
    yticklabels({'M_1','M_2','M_3'})
    ylabel('Model','FontSize',10)
    xlabel('SS periods [-]')
    title('Model choice (Baseline)')

subplot(3,1,2)
    area([0, 50],[5, 5],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2,'LineStyle','none')
    hold on
    area([100, 150],[5, 5],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2,'LineStyle','none')
    area([200, 250],[5, 5],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2,'LineStyle','none')
    
    %%%%%%%%%% 
    % GAMBIS %
    %%%%%%%%%%
    [~,temp] = max(g.probModelArray);
    temp(gThres:end) = zeros(temp3,length(gThres:temp4));
    stairs(0:g.tEnd,temp,'b.','MarkerSize',10)
    
%     text(25,4,'H','HorizontalAlignment','center')
%     text(125,4,'D_1','HorizontalAlignment','center')
%     text(225,4,'D_2','HorizontalAlignment','center')
%     
    ax = gca;
    ax.YGrid = 'on';
    ax.GridLineStyle = '-';
    
    xlim([0, tEnd])
    ylim([0.5, 4.5])
    yticks([1,2,3])
    yticklabels({'M_1','M_2','M_3'})
    ylabel('Model','FontSize',10)
    xlabel('SS periods [-]')
    title('Model with highest probability (ECM)')
    
subplot(3,1,3)
    area([0, 50],[5, 5],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2,'LineStyle','none')
    hold on
    area([100, 150],[5, 5],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2,'LineStyle','none')
    area([200, 250],[5, 5],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2,'LineStyle','none')
    
    %%%%%%%%%% 
    % GAMBIS %
    %%%%%%%%%%
    [~,temp] = max(d.probModelArray);
    temp(dThres:end) = zeros(temp3,length(dThres:temp4));
    stairs(0:d.tEnd,temp,'r.','MarkerSize',10)
    
%     text(25,4,'H','HorizontalAlignment','center')
%     text(125,4,'D_1','HorizontalAlignment','center')
%     text(225,4,'D_2','HorizontalAlignment','center')
%     
    ax = gca;
    ax.YGrid = 'on';
    ax.GridLineStyle = '-';
    
    xlim([0, tEnd])
    ylim([0.5, 4.5])
    yticks([1,2,3])
    yticklabels({'M_1','M_2','M_3'})
    ylabel('Model','FontSize',10)
    xlabel('SS periods [-]')
    title('Model with highest probability (JSDR)')

if pr(1)
    print(f4,'Results_ModelChoice','-dpdf')
end
if pr(2)
    print(f4,'-r1200','-dtiff','Results_ModelChoice.tif');
end

%% 2. Probability vs. Chosen Model B
f5 = figure(5);
leg3 = {'H','D_1','D_2','Shutdown'};

subplot(3,1,1)
    area([0, 50],[5, 5],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2,'LineStyle','none')
    hold on
    area([100, 150],[5, 5],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2,'LineStyle','none')
    area([200, 250],[5, 5],'FaceColor',[0.5,0.5,0.5],'FaceAlpha',0.2,'LineStyle','none')
    
    stairs(1:bl.tEnd,bl.modelArrayProb(2:end),'k.','MarkerSize',10)
    
    text(25,4,'H','HorizontalAlignment','center')
    text(125,4,'D_1','HorizontalAlignment','center')
    text(210,4,'D_2','HorizontalAlignment','center')
    
    ax = gca;
    ax.YGrid = 'on';
    ax.GridLineStyle = '-';
    
    xlim([0, tEnd])
    ylim([0.5, 4.5])
    yticks([1,2,3])
    yticklabels({'M_1','M_2','M_3'})
    ylabel('Model','FontSize',10)
    xlabel('SS periods [-]')
    title('Model choice (Baseline)')

subplot(3,1,2)    
    area(0:g.tEnd,g.probModelArray')
    
    %%%%%%%%%%
    % GAMBIS %
    %%%%%%%%%%
    hold on
    area([gThres temp4],[2 2],'FaceColor','w')

    yticks(0:0.25:1)
    xlim([0, tEnd])
    ylim([0 1])  
   
    legend(leg3,'FontSize',6,'Location','best')
    xlabel('SS periods [-]')
    title('Model probability (ECM)')
    
subplot(3,1,3)
    area(0:d.tEnd,d.probModelArray')
    
    %%%%%%%%%%
    % GAMBIS %
    %%%%%%%%%%
    hold on
    area([dThres temp4],[2 2],'FaceColor','w')
    
    yticks(0:0.25:1)
    xlim([0, tEnd])
    ylim([0 1])  

    %legend(leg3,'FontSize',6,'Location','best')
    xlabel('SS periods [-]')
    title('Model probability (JSDR)')

if pr(1)
    print(f5,'Results_ModelChoice_2','-dpdf')
end
if pr(2)
    print(f4,'-r1200','-dtiff','Results_ModelChoice_2.tif');
end

%% 3. Profit Comparison (together)
f6 = figure(6);
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
    xlim([1,tEnd])
    legend(leg,'Location','best','FontSize',9)

if pr(1)
    print(f6,'Results_Profit','-dpdf')
end
if pr(2)
    print(f6,'-r1200','-dtiff','Results_Profit.tif');
end

%% 3. Profit Comparison (alone)

f7 = figure(7);

%computing the value
    prod_bl = bl.yPlantArray(2,1:3:end)*3600;
    prod_g = g.yPlantArray(2,1:3:end)*3600;
    prod_d = d.yPlantArray(2,1:3:end)*3600;

subplot(3,1,1), 
    plot(1:bl.tEnd + 1,prod_bl,'-kx','MarkerSize',5)
    hold on
    plot(1:g.tEnd + 1,prod_g,'-bo','MarkerSize',5)
    plot(1:d.tEnd + 1,prod_d,'-rd','MarkerSize',5)
%     %hold on 
%     a = area([blThres + 1 dThres + 1],[10000 10000],'FaceColor','y');
%     a.FaceAlpha = 0.1;
%     
    grid on
    title('Total oil production')
    ylabel('Flowrate [m3/h]','FontSize',10)
    xlabel('SS periods [-]','FontSize',10)
    xlim([1,tEnd])
    ylim([0,7500])
    yticks(0:2500:7500)

    leg = {'Baseline','ECM','JSDR'}; 
    legend(leg,'Location','southwest','FontSize',6)
    
%       %creating an inset 
%     % create smaller axes in top right, and plot on it
%     axes('Position',[.2 .5 .3 .15])
%     box on
%     plot(1:g.tEnd + 1,prod_bl,'-kx','MarkerSize',5)
%     hold on
%     plot(1:g.tEnd + 1,prod_g,'-bo','MarkerSize',5)
%     plot(1:d.tEnd + 1,prod_d,'-rd','MarkerSize',5)
%     %hold on 
%     a = area([blThres + 1 dThres + 1],[10000 10000],'FaceColor','y');
%     a.FaceAlpha = 0.1;
%     
%     grid on
%     title('Total oil production (inset)')
%     %ylabel('[m3/h]','FontSize',10)
%     xlabel('SS periods [-]','FontSize',10)
% 
%     ylim([1.45, 1.55])
%     xlim([0 10])
%     

if pr(1)
    print(f7,'Results_Profit_2','-dpdf')
end
if pr(2)
    print(f7,'-r1200','-dtiff','Results_Profit.tif');
end