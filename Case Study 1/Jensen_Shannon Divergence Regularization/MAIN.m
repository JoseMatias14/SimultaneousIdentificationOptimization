% Codes of: Matias et al. Simultaneous Online Model Identification and Control UsingModifier Adaptation and Reinforcement Learning

% Case study 1: Nine different models of the CSTR reactor are available. Difference between the models lies in the reaction kinetic parameter &
% set of reactions (structural mismatch)
%
% Approach used: Jensen_Shannon Divergence Regularization

% Other m-files required: ParametersBlockReactor.m;
% InitialConditionBlockReactor.m; OptimizationBoundsBlockReactor.m;
% PlantModel.m; CDAGradient.m; ECM.m; ReactorModel.m; EconomicOptimizationSS.m; ResultsPlot.m

% MAT-files required: results_standard.mat(from: ...\Baseline Case - Standard Method)
%                     PlantSurface_BR.mat (from:...\Pre Calculation)
%                     covLow.mat           (from:...\Pre Calculation)

clear  
close all
clc

%adding casadi
import casadi.*

% For reproducibility
rng('default')  

%saving file
results = 'results_JSDR';

%% Simulation tuning
%number of ss periods in the simulation
nIter = 50; 

%paramters
par = ParametersBlockReactor;

%initial condition
[dx0,u0] = InitialConditionBlockReactor;

%variable bounds
[~,~,~,ubu] = OptimizationBoundsBlockReactor;

%System parameters for plant model
parPlant = par;

%state output mapping (only the states Xb, Xc & Xd are measured)
H = [0, 1, 0, 0;
     0, 0, 1, 0;
     0, 0, 0, 1];
par.H = H;

%% Online Model Selection + MAy tuning
% number of system measurements
ma.nMeas = 3;
% number of system inputs (only Fa)
ma.nInput = 1; 

%first order filters for the modifiers 
ma.Keps = 0.5*eye(ma.nMeas);        %bias filter (same for all measurements)
ma.Klam = 0.5*eye(ma.nInput);       %gradient filter (same for all measurements) 

%first order input filters  
ma.Kinput = 0.25*eye(ma.nInput);    

% Multimodel tuning
% m1, m2, m3: A + B -k1-> C & 2B -k2-> D & pure D
% m1: k1 = 0.75,  k2 = 1.2 
% m2: k1 = 0.72,  k2 = 1.5 
% m3: k1 = 0.75,  k2 = 1.5  (TRUE MODEL)

% m4, m5, m6: A + B -k3-> 2C & C + 2B -k4-> D & pure D
% m4: k3 = 0.0254, k4 = 0.0166 
% m5: k3 = 0.02,   k4 = 0.02 
% m6: k3 = 0.028,  k4 = 0.02 

% m7, m8, m9: A + B -k1-> C & D as impurity
% m1: k1 = 0.1295
% m2: k1 = 0.135
% m3: k1 = 0.12 

% number of available models
ma.nModels = par.nr; 

%% Gradient Estimation tuning
% Central Difference Approach - CDA
h = 0.75; %perturbation step size

%% Initialization
count = 1; %counter for the number of executions

% Simulation
dxk = dx0;
uk = u0;
%initial guess of the model structure
modelChoicek = 9; 

% uniform prior probability. All the model structure have the same
% a priori probability of being the correct one 
prob_k = 1/par.nr*ones(par.nr,1);

% MAy modifiers are initialized with 0
for j = 1:ma.nModels
    epsk{j} = zeros(ma.nMeas,1); 
    lambdak{j} = zeros(ma.nInput,ma.nMeas); 

end

%% For plotting 
%1. Optimal values
uOptArray = u0;

%2. MAy modifiers 
for j = 1:ma.nModels
    epsArray{j} = epsk{j}; 
    
    for i = 1:ma.nMeas
        lambdaArray{j,i} = lambdak{j}(:,i);
    end  
end

%2. Gradient 
for i = 1:ma.nMeas
    %plant measurements gradients
    gradMeasPlantArray{i} = [];
    
    %plant measurements gradients estimates
    gradPlantHatArray{i} = [];

end

%3. Plant
inputPlantArray = []; % computed by MAy + plant excitation -> contains all the inputs aplied during the simulation
yPlantArray = [];     % plant measurements
yNoNoiseArray = [];   % plant measurements (no noise)

%4. Models
% probability that the model is correct
probModelArray = prob_k;

%which model has the highest probability at the current iteration
modelArrayProb = modelChoicek; 

%store all the inputs used for finding the expected value w.r.t. \lambda
% these inputs are associated with different values of the \lambda
ukLambdaArray = [];

%% Simulation
while count <= nIter
    
     fprintf('     iter. >>> %0.0f \n',count) 
    
    %1. simulate plant (SS) and save plant data
    [~,~,~,~,dxk,yValuePlant,FendGradMeas] = PlantModel(dxk,uk,parPlant);
    
        %%%%%%%%%%%%%%%%% saving plant information
        inputPlantArray = [inputPlantArray, uk];
        yPlantArray = [yPlantArray, H*yValuePlant];
        yNoNoiseArray = [yNoNoiseArray, H*dxk];

        % true gradients
        gradTemp = H*FendGradMeas;
        for i = 1:ma.nMeas
            gradMeasPlantArray{i} = [gradMeasPlantArray{i},gradTemp(i,:)'];
        end
        %%%%%%%%%%%%%%%%%%
        
    %2. estimating gradient using CDA
    [gradPlantHat,uk_h,yk_h,xk_h] = CDAGradient(dxk,uk,parPlant,h,H);
    
        %%%%%%%%%%%%%%%%% saving inputs for plant excitation + system measurements at these points
        inputPlantArray = [inputPlantArray,uk_h];
        yPlantArray = [yPlantArray,yk_h]; 
        yNoNoiseArray = [yNoNoiseArray, xk_h];

        % estimated gradients
        for i = 1:ma.nMeas
            gradPlantHatArray{i} = [gradPlantHatArray{i}, gradPlantHat(i,:)'];
        end
        
         %indicating that this is a probing trial for gradient estimation
         modelArrayProb = [modelArrayProb,0,0]; 
        %%%%%%%%%%%%%%%%%  

    %3. Optimizing 
    [pi_k,prob_k,uk_lambda,epsk,lambdak] = JSDR(dxk,uk,par,yValuePlant,gradPlantHat,epsk,lambdak,prob_k,H,ma);    

        % implementing input filter (for all lambda's in the expectation)
        uk_lambda = uk + ma.Kinput*(uk_lambda - uk); 
        ukLambdaArray = [ukLambdaArray, uk_lambda];
       
        % actual implemented value
        uk = uk + ma.Kinput*(pi_k - uk); 

        % probability (a posteriori | after the measurements) associated with each of the ma.nModels models
        probModelArray = [probModelArray, prob_k];
        
        % analyzing which model has the highest probability of being the
        % correct one
        [~,modelChoicek] = max(prob_k);        
        
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
[~,~,~,~,dxk,yValuePlant,FendGradMeas] = PlantModel(dxk,uk,parPlant);

    %%%%%%%%%%%%%%%%% plant information
    inputPlantArray = [inputPlantArray, uk];
    yPlantArray = [yPlantArray, H*yValuePlant];
    yNoNoiseArray = [yNoNoiseArray, H*dxk];

%% Saving the results and ploting
save(results,'uOptArray','epsArray','lambdaArray','gradPlantHatArray','gradMeasPlantArray','yPlantArray','yNoNoiseArray','inputPlantArray','ukLambdaArray','probModelArray','modelArrayProb','nIter','par','ma','parPlant');

%plotting the results
ResultsPlot