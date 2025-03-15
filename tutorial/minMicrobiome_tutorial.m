% Tutorial on minimal microbiome analysis 
% Install the minMicrobiome.m and createCommModel.m from the main GitHub page and add the path to these files

% Initialize diet
%dietFilePath = '/path/to/diet/file';
dietConstraints=readtable(dietFilePath);

%% (Optional) To modify the diet constraints to match the model 
dietConstraints=table2cell(dietConstraints);
dietConstraints(:,2)= cellstr(num2str(cell2mat(dietConstraints(:,2))));
dietCom = dietConstraints;
dietCom (:,1) = strrep(dietCom(:,1),'(e)','(u)'); 

%% Creating a community model
modelList = table2cell(readtable('/path/to/modelListfile'));

%modelPath = '/path/to/AGORA/matfiles'; 
modelCell = {}; spBm = {}; spATPM = {}; options = struct();
for k = 1: length(modelList)
    model = readCbModel([modelPath filesep modelList{k,1} '.mat']);
    modelCell(k) = {model};
    spBm(k,1) = {model.rxns{strncmp(model.rxns,'biomass',7)}};
    spATPM(k,1) = {'DM_atp_c_'};
end
options.spBm = spBm;
options.spATPM = spATPM;
options.sepUtEx = false;

% Building community model
[modelCom, infoCom, indCom, msg] = createCommModel(modelCell, options);

% Applying diet to the community
changeCobraSolver('ibm_cplex', 'LP');
modelCom=useDiet(modelCom,dietCom);
modelCom.lb(startsWith(modelCom.rxns,'EX_')& ...
    contains(modelCom.rxns,'_org')) = -1000;

%% Parameters for minimal microbiome

changeCobraSolver('ibm_cplex', 'MILP');
solutionPert = {};
options.gr_frac = 0.8;
options.scfa_frac = 0.8;
options.gr_opt_frac = 0.99; % Fraction of growth rate for FVA on exchange reactions
options.met_names = {'EX_for(u)';'EX_but(u)'; 'EX_ppa(u)'}; %format, butyrate and propionate 'EX_ac(u)'; 
options.constraint = 1; %1: sum SCFA, 2: acetate, 3: butyrate, 4: propionate
options.maxMILP = 8; %Provide maximum of 10 for resonable computation time
options.iter = 50; %Trade-off between maxMILP and iter

tic
[solutionPert, supp_out] = MinMicrobiome(modelCom, options);
toc

% Another example on changing the optimal fraction of growth rate
options2 = options;
options2.gr_opt_frac = 0.9;
tic
[solutionPert_2, supp_out_2] = MinMicrobiome(modelCom, options2);
toc
