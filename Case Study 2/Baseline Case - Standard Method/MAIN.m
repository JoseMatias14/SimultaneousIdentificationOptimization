% Codes of: Matias et al. Simultaneous Online Model Identification and Control UsingModifier Adaptation and Reinforcement Learning

% Case study 2: three different models (healthy/degraded 1 and degraded 2)
% of a compressor in a subsea compression station are available. Difference
% between the models lies in the the function for computing the compressor
% efficiency
%
% Approach used: Baseline 

% Other m-files required: ParametersBlockReactor.m;
% InitialConditionBlockReactor.m; OptimizationBoundsBlockReactor.m;
% PlantModel.m; CDAGradient.m; ECM.m; ReactorModel.m; EconomicOptimizationSS.m; ResultsPlot.m

% MAT-files required: results_standard.mat(from: ...\Baseline Case - Standard Method)
%                     PlantSurface_BR.mat (from:...\Pre Calculation)
%                     covMatrix.mat           (from:...\Pre Calculation)

clear
close all
clc

%adding casadi
import casadi.*

% For reproducibility
rng('default')  

%saving file
results = 'results_standard';


%% Simulation tuning
%number of ss periods in the simulation
tEnd = 250;

%paramters
par = ParametersSubseaGas;

%initial condition
[x0,u0,f0] = InitialConditionSubseaGas(par);

%variable bounds
[lbx,ubx,lbu,ubu] = OptimizationBoundsSubseaGas;

%System parameters for plant model
parPlant = par;

%state output mapping
H = zeros(8,10);
H(1,1) = 1; %choke outlet pressure: P_choke_out*1e7 [Pa] -- NOTE THAT the values are scaled!
H(2,2) = 1; %volumetric flow from the reservoir: q_res [m3/s]
H(3,3) = 1; %inlet compressor volumetric flowrate: q_comp_in [m3/s]
H(4,4) = 1; %compressor outlet temperature: T_comp_out*1e2 [K]
H(5,5) = 1; %compressor outlet pressure: P_comp_out*1e7 [Pa]
H(6,8) = 1; %compressor power: Pow*1e4 [W]
H(7,9) = 1;%inlet pump volumetric flowrate: q_pump_in*1e-1 [m3/s]
H(8,10) = 1; %pump power: Pop*1e3 [W]
par.H = H;

%% Online Model Selection + MAy tuning
% number of system measurements
ma.nMeas = 8;

% number of system inputs (choke opening and compressor speed)
ma.nInput = 2;

%first order filters for the modifiers
ma.Keps = 0.75*eye(ma.nMeas);   %bias filter (same for all measurements)
ma.Klam = 0.5*eye(ma.nInput);	%gradient filter (same for all measurements) 

%first order input filters  
ma.Kinput = 0.25*eye(ma.nInput); 

% Multimodel tuning
% m1: healthy compressor
% m2: degraded compressor
% m3: heavily degraded compressor

% number of available models
ma.nModels = 3; 

%% Gradient Estimation tuning
% Central Difference Approach - CDA
h = 0.05;%perturbation size

%% Initialization
% simulation counter
count = 1;

% since the plant behavior changes with time (healthy->degraded 1-> degraded 2), we created a counter for the plant states. It is re.initialized everytime the compressor is repaired
countPlant = count; 

% Simulation
xk = x0;
fk = f0; % note that in this case study, we have a feed stream
uk = u0;

%initial guess of the model structure
modelChoicek = 1; %healthy model 

% uniform prior probability. All the model structure have the same
% a priori probability of being the correct one 
prob_k = 1/ma.nModels*ones(ma.nModels,1);

%MA - modifiers are considered as if the first point has no mismatch
for j = 1:ma.nModels
    epsk{j} = zeros(ma.nMeas,1);
    lambdak{j} = zeros(ma.nInput,ma.nMeas);

end

%the feed pressure oscillates. It is approximated by a square wave. Here we
%create an auxiliary variable with the shape of oscillation 
% 0.1: adjusts the oscillation frequency 
tAux = 1:tEnd;
sqAux = square(0.1*tAux);

%% For plotting - several values are stored for ploting
%1. Optimal values
uOptArray = u0;

%2. MAy modifiers 
for j = 1:ma.nModels
    epsArray{j} = epsk{j}; 

    for i = 1:ma.nMeas
        lambdaArray{j,i} = lambdak{j}(:,i);
    end
end

%3. Gradient
for i = 1:ma.nMeas
    %plant measurements gradients
    gradMeasPlantArray{i} = [];
    
    %plant measurements gradients estimates
    gradPlantHatArray{i} = [];

end

%4. Plant
inputPlantArray = []; % computed by MAy + plant excitation -> contains all the inputs aplied during the simulation
yPlantArray = [];     % plant measurements

%5. Models
% inputs computed using the available models
uModelArray = [];

% probability that the model is correct
probModelArray = prob_k;

%which model has the highest probability at the current iteration
modelArrayProb = modelChoicek; 

%% Simulation
count = 1;

while count <= tEnd

    fprintf('     iter. >>> %0.0f \n',count)
    %1. simulate plant (SS) and store plant data
    
    % updating the feed stream pressure (oscillates according to square
    % wave
    if sqAux(count) == 1
        fk(3) = 1.1;
    else %sqAux == -1
        fk(3) = 1.2;
    end
    
    % simulating the plant model
    [xk,xmk,gradPlant] = PlantModel(xk,uk,fk,parPlant,count);

        %%%%%%%%%%%%%%%%% saving plant information
        inputPlantArray = [inputPlantArray, uk];
        yPlantArray = [yPlantArray, H*xmk];

        % true gradients - we don't save all gradients, only the ones
        % related to the measured values
        gradTemp = H*gradPlant;
        for i = 1:ma.nMeas
            gradMeasPlantArray{i} = [gradMeasPlantArray{i},gradTemp(i,:)'];
        end
        %%%%%%%%%%%%%%%%%%

    %2. Gradient estimation using CDA
    [gradPlantHat,uk_h1,uk_h2,yk_h1,yk_h2] = CDAGradient(H*xmk,xk,uk,fk,parPlant,h,H,count);

        %%%%%%%%%%%%%%%%% saving inputs for plant excitation + system measurements at these points
        inputPlantArray = [inputPlantArray,uk_h1,uk_h2];
        yPlantArray = [yPlantArray,yk_h1(1:ma.nMeas),yk_h2(1:ma.nMeas)];

        for i = 1:ma.nMeas
            gradPlantHatArray{i} = [gradPlantHatArray{i},gradPlantHat(i,:)'];
        end
        %%%%%%%%%%%%%%%%

   %3. Optimizing 
    [pi_k,~,~,pi_k_index,epsk,lambdak] = OnlineMaintenance(xk,uk,fk,par,H*xmk,gradPlantHat,epsk,lambdak,prob_k,H,ma);    
        
        %choosing the model with the highest probability
        modelChoicek = pi_k_index;

   %4. Updating the plant inputs
    %choosing the optimal input compute when model == modelChoice 
    uStar = pi_k;

    % implementing input filter
    uk = uk + ma.Kinput*(uStar - uk);
        
   %5. For plotting
        % Optimization
        uOptArray = [uOptArray, uk];

        % MAy
        for j = 1:ma.nModels
            epsArray{j} = [epsArray{j}, epsk{j}];

            for i = 1:ma.nMeas
                lambdaArray{j,i} = [lambdaArray{j,i}, lambdak{j}(:,i)];
            end
        end

        % Model choice 
        modelArrayProb = [modelArrayProb,modelChoicek];

    %loop
    count = count + 1;

end

%calculating the results for the last point
% simulate plant (SS) and save plant data
[xk,xmk,gradPlant] = PlantModel(xk,uk,fk,parPlant,count);

    %%%%%%%%%%%%%%%%% plant information
    inputPlantArray = [inputPlantArray, uk];
    yPlantArray = [yPlantArray, H*xmk];

%% Saving the results and ploting
save(results,'inputPlantArray','yPlantArray','gradMeasPlantArray','gradPlantHatArray',...
              'uOptArray','uModelArray','probModelArray','uChosenModelArray','modelArrayProb','modelArrayU',...
              'tEnd','ma','par','parPlant');

ResultsPlotting          
          
beep