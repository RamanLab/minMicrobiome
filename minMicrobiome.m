function [solutionPert, supp_out] =  minMicrobiome(modelCom,options)
milp_max = 50;
tol = 1E-3;
n_org_actual = size(modelCom.modelID,1); %No. of organisms in original community
%% The optional inputs assignment
if isfield(options,'gr_frac')
    gr_frac = options.gr_frac;
else
    gr_frac = 0.8;
end
if isfield(options,'scfa_frac')
    scfa_frac = options.scfa_frac;
else
    scfa_frac = 0.8;
end

if isfield(options,'c')
    c = options.c;
else
    c = 1000;
end
if isfield(options,'u')
    u = options.u;
else
    u = 0.01;
end

if isfield(options, 'constraint')
    cnstrt = options.constraint;
else
    cnstrt = 1;
end

if isfield(options, 'maxMILP')
    maxMILP = options.maxMILP;
else
    maxMILP = 10;
end

if isfield(options, 'iter')
    iter = options.iter;
else
    iter = 3;
end

if isfield(options, 'max_del_rounds')
    max_del_rounds = options.max_del_rounds;
else
    max_del_rounds = 10;
end

if isfield(options, 'gr_opt_frac')
    gr_opt_frac = options.gr_opt_frac;
else
    gr_opt_frac = 0.99;
end

if (isfield(options, 'constraint_wt') && cnstrt == 1)
    wtg = options.constraint_wt;
else
    wtg = 1;
end

if isfield(options, 'met_names') 
    scfa_pat = options.met_names;
else
    scfa_pat = {'EX_ac(u)'; 'EX_but(u)'; 'EX_ppa(u)'};  %acetate, butyrate and propionate
end

solutionPert = cell(1,iter);
supp_out = struct();
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
disp(max_scfa);
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
    supp_out.del_seq(i_val,:) = del_seq;
    
    
    while(l>maxMILP && ii<=n_org_actual && del_rounds<max_del_rounds)
        if (all(del_list~=del_seq(ii)))
            pat = strcat('org',int2str(del_seq(ii)));
            n_r = find(cellfun(@(x) (length(char(x))>length(pat)) ...
                && strcmpi(pat,x(length(char(x))-(length(pat)-1):end)),modelCom.rxns));
            modelCom_new = changeRxnBounds(modelCom_new, modelCom_new.rxns(n_r),0, 'b');
            
            res_new = optimizeCbModel(modelCom_new);
            if (res_new.f>= gr_frac* gr_max)
                pertModel_scfa = modelCom_new;
                for k = 1:length(bm)
                    pertModel_scfa = addCOBRAConstraints(pertModel_scfa, ...
                        bm_rxns(k), gr_opt_frac*res_new.v(bm(k)),'dsense', 'G');
                end
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
                end
            else
                modelCom_new = pertModel;
            end
            ii = ii+1;
        else
            ii = ii+1;
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
        disp('No organism was removed')
        res_new = res_Com;
    else
        org_thrown_out{i_val} = modelCom.modelID(del_list);
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
    MILPproblem.x0 = [res_new.v; ones(n_org_actual,1)];
    
    solutionPert{i_val} = solveCobraMILP(MILPproblem);
    
    if exist('org_thrown_out')
        supp_out.org_thrown_out = org_thrown_out;
    end
    
end
if all(cellfun(@isempty,solutionPert))
    disp('Try again');
end
