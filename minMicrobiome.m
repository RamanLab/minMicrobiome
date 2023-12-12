function [solutionPert, supp_out] =  minMicrobiome(modelCom,options)
startTime = tic;    %Starting the timer
milp_max = 100;     %Max possible size for MILP - change according to requirement


%Initializations
tol = 1E-3;         %Used to determine SCFA production rates
n_org_actual = size(modelCom.modelID,1); %No. of organisms in original community
supp_out = struct();
solutionPert = {};
supp_out.org_thrown_out = {};
supp_out.num_min_orgs = [];
supp_out.del_seq = [];
supp_out.scfa_minMicrobiome = [];
supp_out.minMicrobiomes = {};
orgs_in_minMicrobiome = {};
exitOuterLoop = false;

%% The optional inputs assignment
if isfield(options,'gr_frac')   %Fraction of growth rate needed to be preserved in the minimal microbiome
    gr_frac = options.gr_frac;
else
    gr_frac = 0.8;
end
if isfield(options,'scfa_frac') %Fraction of SCFA needed to be preserved in the minimal microbiome
    scfa_frac = options.scfa_frac;
else
    scfa_frac = 0.8;
end

if isfield(options,'c') %For coupling constraints to avoid fluxes without biomass
    c = options.c;
else
    c = 1000;
end
if isfield(options,'u') %For coupling constraints to avoid fluxes without biomass
    u = options.u;
else
    u = 0.01;
end

if isfield(options, 'constraint')   %The functionality constraint
    cnstrt = options.constraint;
else
    cnstrt = 1;
end

if isfield(options, 'maxMILP')      %No. of members (integer variables) in MILP that is to be minimized
    maxMILP = options.maxMILP;
else
    maxMILP = 10;
end

if isfield(options, 'MILP_runs')    %No. of times MILP will be run by changing the initial conditions
    MILP_runs = options.MILP_runs;
    if (MILP_runs > n_org_actual)
        MILP_runs = 1;
        fprintf("Since MILP_runs > max no. of organisms, it has been set to %d.\n", n_org_actual);
    end
else
    MILP_runs = n_org_actual;
end

if isfield(options, 'iter')     %No. of iterations -- staring from the deletion step
    iter = options.iter;
else
    iter = 3;
end

if isfield(options, 'time_budget')  %How much time can we afford to run?
    time_budget = options.time_budget;
else
    time_budget = 10800;
end

if isfield(options, 'max_del_rounds')   %No. of times deletion sequence is repeated to arrive at 'maxMILP' elements in the MILP
    max_del_rounds = options.max_del_rounds;
else
    max_del_rounds = 10;
end

if isfield(options, 'gr_opt_frac')  %Optimum fraction of biomass growth for SCFA production
    gr_opt_frac = options.gr_opt_frac;
else
    gr_opt_frac = 0.99;
end

if (isfield(options, 'constraint_wt') && cnstrt == 1)   %Weightage on functionality constraint
    wtg = options.constraint_wt;
else
    wtg = 1;
end

if isfield(options, 'met_calc_minMicrobiome')  %Should we calculate the metabolite rate from the resultant minimal microbiome?
    met_calc_min = options.met_calc_minMicrobiome;
else
    met_calc_min = 'yes';
end

if isfield(options, 'met_names')
    scfa_pat = options.met_names;   %Names of the metabolites of interest - to be used in the functionality
else
    scfa_pat = {'EX_ac(u)'; 'EX_but(u)'; 'EX_ppa(u)'};  %acetate, butyrate and propionate
end

%% Coupling constraints
for i = 1:n_org_actual
    pat = strcat('org',int2str(i));
    n_r = find(cellfun(@(x) (length(char(x))>length(pat)) ...
        && strcmpi(pat,x(length(char(x))-(length(pat)-1):end)),modelCom.rxns));
    idx = find(strncmp(modelCom.rxns(n_r), 'biomass',7));
    n_rB = n_r(idx);
    n_r(idx) = [];
    modelCom = coupleRxnList2Rxn(modelCom,modelCom.rxns(n_r), ...
        modelCom.rxns(n_rB), c, u);
end
res_Com = optimizeCbModel(modelCom);

%% SCFA and biomass constraints
switch cnstrt
    case 1
        for k=1:length(scfa_pat)
            scfa_list(k) = find(strcmp(modelCom.rxns, scfa_pat{k}));
        end
    case 2
        scfa_list = find(strcmp(modelCom.rxns, scfa_pat{1}));
    case 3
        scfa_list = find(strcmp(modelCom.rxns, scfa_pat{2}));
    case 4
        scfa_list = find(strcmp(modelCom.rxns, scfa_pat{3}));
end
bm = find(modelCom.c);
bm_rxns = modelCom.rxns(bm);
gr_max = res_Com.f; gr_ind = res_Com.v(bm);

%% Choose the constraint
%Find max SCFAs possible from the community at known growth rate
pertModel_scfa = modelCom;
for k = 1:length(bm)
    pertModel_scfa = addCOBRAConstraints(pertModel_scfa, ...
        bm_rxns(k), gr_opt_frac* gr_ind(k),'dsense', 'G');
end

pertModel_scfa.c(:) = 0;
pertModel_scfa.c(scfa_list) = wtg;
LPproblem = buildLPproblemFromModel(pertModel_scfa);
solutionPert_scfa = solveCobraLP(LPproblem);
max_scfa = solutionPert_scfa.full(scfa_list);
fprintf('Maximum SCFA production by the large community is %f\n',max_scfa);
switch cnstrt
    case 1
        if sum(max_scfa)<=tol
            fprintf('Net SCFA production is 0! \nConstraint added when SCFA exchange is positive.\n')
        end
    case 2
        if max_scfa<=tol
            disp('Community does not produce acetate');
            return
        end
    case 3
        if max_scfa<=tol
            disp('Community does not produce butyrate');
            return
        end
    case 4
        if max_scfa<=tol
            disp('Community does not produce propionate');
            return
        end
end

supp_out.gr_max = gr_max;
supp_out.gr_ind = gr_ind;
supp_out.max_scfa = max_scfa;
for i_val = 1:iter
    if isfield(options,'del_seq')
        del_seq = options.del_seq(i_val,:);
    else
        del_seq = randperm(n_org_actual,n_org_actual);
    end
    del_rounds = 1;
    del_list = [];
    l = n_org_actual; ii =1;
    modelCom_new = modelCom;
    pertModel = modelCom;
    
    while(l>maxMILP && ii<=n_org_actual && del_rounds<max_del_rounds)
        if (all(del_list~=del_seq(ii)))
            pat = strcat('org',int2str(del_seq(ii)));
            n_r = find(cellfun(@(x) (length(char(x))>length(pat)) ...
                && strcmpi(pat,x(length(char(x))-(length(pat)-1):end)),modelCom.rxns));
            modelCom_new = changeRxnBounds(modelCom_new, modelCom_new.rxns(n_r),0, 'b');
            
            res_new = optimizeCbModel(modelCom_new);
            if (res_new.f>= gr_frac* gr_max)
                pertModel_scfa = modelCom_new;
                %Growth rate constraint - individual gr rates considered
                for k = 1:length(bm)
                    pertModel_scfa = addCOBRAConstraints(pertModel_scfa, ...
                        bm_rxns(k), gr_opt_frac*res_new.v(bm(k)),'dsense', 'G');
                end
                % Growth rate constraint based on overall gr rate. Less strict
                % than the previous one - Final result could contain more
                % number of species in minMicrobiome when this constraint
                % is used.
                %                 pertModel_scfa = addCOBRAConstraints(pertModel_scfa, ...
                %                           bm_rxns, gr_frac*gr_max, 'dsense', 'G');
                
                pertModel_scfa.c(:) = 0;
                pertModel_scfa.c(scfa_list) = wtg;
                LPproblem = buildLPproblemFromModel(pertModel_scfa);
                solutionPert_scfa = solveCobraLP(LPproblem);
                if solutionPert_scfa.stat ==1
                    scfa_new = solutionPert_scfa.full(scfa_list);
                else
                    ii = ii+1;
                    continue;
                end
                
                if all((max_scfa>0 & scfa_new>=scfa_frac*max_scfa)|(max_scfa<=0))
                    pertModel = modelCom_new;
                    l = l-1;
                    del_list(n_org_actual-l) = del_seq(ii);
                    
                else
                    modelCom_new = pertModel;
                    fprintf('Organism %d cannot be deleted! Trying another one...\n',del_seq(ii));
                    
                end
            else
                modelCom_new = pertModel;
                fprintf('Organism %d cannot be deleted! Trying another one...\n',del_seq(ii));
                
            end
            ii = ii+1;
            if(length(del_seq)<ii)
                break
            end
        else
            ii = ii+1;
            if(length(del_seq)<ii)
                break
            end
        end
        
        
        if (l>milp_max && ii>n_org_actual)
            ii = 1;
            del_rounds = del_rounds+1;
        end
    end
    if l>milp_max
        fprintf('MILP larger than %d!\n', milp_max);
        continue
    end
    if isempty(del_list)
        fprintf('No organism was removed.\n')
        res_new = res_Com;
    else
        org_thrown_out = modelCom.modelID(del_list);
    end
    %% Add new variables & constraints
    %Binary variables
    new_var = strcat('binary_', modelCom.modelID);
    pertModel = addCOBRAVariables(pertModel,new_var,'lb',0,'ub',1);
    for i = 1:n_org_actual
        pat = strcat('org',int2str(i));
        n_r = find(cellfun(@(x) (length(char(x))>length(pat)) ...
            && strcmpi(pat,x(length(char(x))-(length(pat)-1):end)),modelCom.rxns));
        idx = find(strncmp(modelCom.rxns(n_r), 'biomass',7));
        n_rB = n_r(idx);
        
        pertModel = addCOBRAConstraints(pertModel, {modelCom.rxns{n_rB}, ...
            new_var{i}}, 0,'c', [-1, modelCom.lb(n_rB)], 'dsense', 'L');
        pertModel = addCOBRAConstraints(pertModel, {modelCom.rxns{n_rB}, ...
            new_var{i}}, 0,'c', [-1, modelCom.ub(n_rB)], 'dsense', 'G');
        
        if any(del_list==i)
            pertModel = addCOBRAConstraints(pertModel, new_var{i}, 0, 'dsense', 'E');
        end
    end
    
    % SCFA constraint
    for i = 1:length(scfa_list)
        if max_scfa(i)>0
            scf_c = max_scfa(i) * scfa_frac;
            pertModel = addCOBRAConstraints(pertModel,modelCom.rxns(scfa_list(i)), ...
                scf_c, 'dsense', 'G');
        end
    end
    
    gr_min = gr_frac* res_Com.f;
    pertModel = addCOBRAConstraints(pertModel, bm_rxns, gr_min, 'dsense', 'G');
    
    %% Change objective function - minimize membership vector
    pertModel.c(:) = 0;
    pertModel.evarc(cellfun(@length,regexp(pertModel.evars,'^binary_'))==1) = 1;
    pertModel.osenseStr = 'min';
    
    %% MILP solution
    MILPproblem = buildLPproblemFromModel(pertModel);
    MILPproblem.vartype = char(ones(1,size(modelCom.S,2))*'C');
    MILPproblem.vartype(size(modelCom.S,2)+1: size(modelCom.S,2)+n_org_actual) = 'B';
    
    %Multiple runs changing the starting values of MILP
    for kk = 1:MILP_runs
        int_var_arr = ones(n_org_actual,1);
        if ismember(kk, del_list)
            %if (kk ~= MILP_runs || ~all(cellfun(@isempty,solutionPert)))
            if (iter >1 || (kk~=1 && iter==1))
                fprintf("Run %d of iteration %d: Organism %d has already been deleted. The solution from this run will be redundant. Proceeding to the next run.\n", kk, iter, kk);
                continue;
            end
        end
        int_var_arr(1:kk) = 0;
        MILPproblem.x0 = [res_new.v; int_var_arr];
        solutionPert_in_fun = {};
        solutionPert_in_fun = solveCobraMILP(MILPproblem);
        
        if (solutionPert_in_fun.stat ~=1)
            fprintf('Model %d in iteration %d is infeasible\n', kk, i_val);
            continue;
        end
        
        %Populating solutionPert - redundant solutions are avoided
        repeat_solution = 0;
        if exist ('solutionPert','var')
            for k = 1:size(solutionPert,2)
                if (solutionPert_in_fun.int == solutionPert{k}.int)
                    repeat_solution = 1;
                end
            end
        end
        
        if(repeat_solution == 0)
            solutionPert{end+1} = solutionPert_in_fun;
            if exist('org_thrown_out', 'var')
                supp_out.org_thrown_out(end+1,:) = org_thrown_out;
            end
            if ~isempty(del_list)
                supp_out.del_seq(end+1,:) = del_seq;
            end
            supp_out.num_min_orgs(end+1) = nnz(solutionPert_in_fun.int);
            
            orgs_in_minMicrobiome = modelCom.modelID(find(solutionPert_in_fun.int));
            orgs_stitched = strjoin(orgs_in_minMicrobiome,',\n');
            supp_out.minMicrobiomes{end+1} = orgs_stitched;
        else
            fprintf('Solution %d from iteration %d is redundant\n',kk, i_val);
            continue;
        end
        
        
        
        %% Butyrate production rate of the minimal community
        if strcmp(met_calc_min,'yes')
            model_scfa_calc = modelCom;
            for i = 1:n_org_actual
                pat = strcat('org',int2str(i));
                n_r = find(cellfun(@(x) (length(char(x))>length(pat)) ...
                    && strcmpi(pat,x(length(char(x))-(length(pat)-1):end)),model_scfa_calc.rxns));
                idx = find(strncmp(model_scfa_calc.rxns(n_r), 'biomass',7));
                n_rB = n_r(idx);
                
                if (solutionPert_in_fun.int(i) == 0)
                    
                    for jj = 1:length(n_r)
                        model_scfa_calc = addCOBRAConstraints(model_scfa_calc, n_r(jj), 0, 'dsense', 'E');
                    end
                    
                end
            end
            gr_min = gr_frac* res_Com.f;
            model_scfa_calc = addCOBRAConstraints(model_scfa_calc, bm_rxns, gr_min, 'dsense', 'G');
            
            model_scfa_calc.c(:) = 0;
            model_scfa_calc.c(scfa_list) = wtg;
            LPproblem = buildLPproblemFromModel(model_scfa_calc);
            solutionPert_model_scfa_calc = solveCobraLP(LPproblem);
            
            supp_out.scfa_minMicrobiome(end+1) = solutionPert_model_scfa_calc.obj;
            
            elapsedTime = toc(startTime);
            if (elapsedTime > time_budget)
                disp('Time budget exceeded! Stopping execution! Results based on the run so far are displayed');
                exitOuterLoop = true;
                supp_out.elapsedTime = elapsedTime;
                break;
            end
        else
            elapsedTime = toc(startTime);
            if (elapsedTime > time_budget)
                disp('Time budget exceeded! Stopping execution! Results based on the run so far are displayed.');
                exitOuterLoop = true;
                supp_out.elapsedTime = elapsedTime;
                break;
            end
            
            continue;
        end
        
    end
    if (exitOuterLoop == true)
        break
    end
    
end

%% Wrapping up
if all(cellfun(@isempty,solutionPert))
    disp('Try again');
    elapsedTime = toc(startTime);
    fprintf('Elapsed time is %g seconds\n',elapsedTime);
else
    [supp_out.smallest_minMicrobiome, idx] = min(supp_out.num_min_orgs);
    fprintf('Here is a microbiome with the minimal organisms:\n');
    disp(modelCom.modelID(find(solutionPert{idx}.int)));
    if strcmp(met_calc_min,'yes')
        fprintf('The max butyrate production by this community is %f. \n', supp_out.scfa_minMicrobiome(idx));
        [supp_out.best_scfa_producer_minMicrobiome, idx] = max(supp_out.scfa_minMicrobiome);
        fprintf('The best butyrate producer minimal microbiome is:\n');
        disp(modelCom.modelID(find(solutionPert{idx}.int)))
        fprintf('The max butyrate production by this community is %f. \n', supp_out.scfa_minMicrobiome(idx));
    end
    supp_out.pert_gr = sum(solutionPert{1}.cont(bm));
    supp_out.pert_scfa = solutionPert{1}.cont(scfa_list);
    
    elapsedTime = toc(startTime);
    supp_out.elapsedTime = elapsedTime;
    fprintf('Elapsed time is %g seconds\n',elapsedTime);
end
