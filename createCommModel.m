function [modelCom,infoCom,indCom, msg] = createCommModel(modelCell, options)
%Create a community COBRA model. The model has an extra compartment [u] 
%for inter-organism and community exchange, i.e.,
% (environment) <=> [u] <=> [e_organism1]
%                       <=> [e_organism2] ...
%
%[modelCom,infoCom,indCom] = createCommModel(modelCell,options)
%
%INPUT
% modelCell:    Cell array of COBRA model (e.g., {model1, model2, model3})
%               or a structure with modelCell.org being the model for
%               organism org
% (if a model in 'modelCell' has the field 'metComs', the name in
% model.metComs would be used to map community metabolites instead of the
% original name in model.mets [recommended to provide])
%
% options:      Structure containing the following fields
%     spBm:     Cell array of reaction names of the biomass reaction in
%               each model in 'modelCell'
%   (below are optional, but recommended to provide)
%     spAbbr:   Cell array of abbrivation for each organism in modelCell
%     spName:   Full names of the species
%     spATPM:   Cell array of reaction names of the ATP maintenance reaction in
%               each model in 'modelCell'
%     ('sp' originally for 'species')
%     sepUtEx:  true to have separate uptake and export reactions for community exchange
%               false to have one single exchange reaction for each community metabolite (default true)
%     addExRxns: Add exchange reactions for extracellular metabolites having 
%                no exchange reactions in the model. Otherwise those
%                metabolites will not be connected to the community space (default true)
%     metExId:  Identifier for extracellular metabolites metabolites (default '[e]')
%               If input is the empty string (''), find all metabolites that have exchange reactions
%     rxnField: Cell array of field names in the models that have the same
%               size as rxns. Default to include all fields starting
%               with 'rxn' and include 'grRules', 'rules', 'confidenceScores', 'subSystems'
%     metField: Cell array of field names in the models that have the same size as mets.
%                Default to include all fields starting with 'met' and the field 'csense'
%     geneField: Cell array of field names in the models that have the same size as genes. 
%               Default to include all fields starting with 'gene' and the field 'proteins'
%     compField: Cell array of field names in the models that have the same size as comps (compartments).
%               Default to include all fields starting with 'comp'
%     verbFlag: true to print messages and warnings (default true)
%
%OUTPUT
% modelCom:     COBRA community model with the following extra fields
%               'infoCom' and 'indCom', which are also outputed as separate
%               output argument. See below.
% infoCom:      A structure of useful reaction names and organism names:
%   spBm:          spBm{k} the biomass reaction of info.spAbbr{k}                                        
%   spATPM:        spATPM{k} the ATP maintenance reaction of info.spAbbr{k}                                            
%   rxnSD:         all sink and demand reactions other than the community exchange reactions (EXcom)
%   EXcom:         EXcom{i,1} the community uptake reaction for community metabolite info.metCom{i}
%                  EXcom{i,2} the community production reaction for community metabolite info.metCom{i}
%   EXsp:          EXsp{i,k} the species-community exchange reaction for community metabolite info.metCom{i}
%   Mcom:          community metabolites
%   Msp:           Msp{i,k} the extracellular met of organism k
%                  corresponding to community metabolite Mcom{i}
%   spAbbr:        species abbreviation, used in mets and rxns                       
%                  to identify species-specific metabolites                          
%   spName:        full name of each species        
%   rxnSps:        rxnSps{j} is the abbreviation of species k (spAbbr{k}) 
%                  if reaction j is of species k. 'com' for community exchange reactions    
%   metSps:        metSps{i} is the abbreviation of species k (spAbbr{k}) 
%                  if metabolite i is of species k. 'com' for community exchange reactions 
% indCom:       The corresponding reaction IDs and organism IDs:
%   spBm:          reaction IDs for info.spBm
%   spATPM:        reaction IDs for info.spATPM
%   rxnSD:         reaction IDs for info.rxnSD
%   EXcom:         reaction IDs for info.EXcom                               
%   EXsp:          reaction IDs for info.EXsp. EXsp(i,k) = 0 if info.EXsp{i,k} is empty          
%   Mcom:          metabolite IDs for info.Mcom
%   Msp:           metabolite IDs for info.Msp
%   rxnSps:        index.rxnSps(j) = k implies that reaction j is of                      
%                  species info.spName{k}. 0 for community exchange reactions    
%   metSps:        index.metSps(i) = k implies that metabolite i is of                    
%                  species info.spName{k}. 0 for community metabolites             
% msg:          Messages printed during the model compilation

%% arguement checking
if ~exist('options', 'var')
    options = struct();
end
%get parameters
[spAbbr,spName,spBm,spATPM,metExId,rxnField,metField,geneField,compField,sepUtEx,addExRxns,verbFlag] = getCobraComParams(...
    {'spAbbr','spName','spBm','spATPM','metExId','rxnField','metField','geneField','compField','sepUtEx','addExRxns','verbFlag'}, options);
%organisms' abbreviations and names
nameSpecies = false;
if isstruct(modelCell)
    if isempty(spAbbr)
        spAbbr = fieldnames(modelCell);
    end
    modelCell = struct2cell(modelCell);
else
    if isempty(spAbbr)
        nameSpecies = true;
    end
end
nSp = numel(modelCell);
if nameSpecies
    spAbbr = strcat('org',strtrim(cellstr(num2str((1:nSp)'))));
end
if isempty(spName)
    spName = spAbbr;
end
%biomass reactions
if isempty(spBm)
    error('Please provide the names of the biomass reactions in options.spBm'); 
elseif numel(spBm) ~= nSp
    error('Number of entries in options.spBm not equal to the number of models.');
elseif iscell(spBm)
    rxnBiomassID = zeros(nSp, 1);
    for j = 1: nSp
        rxnBiomassID(j) = findRxnIDs(modelCell{j}, spBm{j});
    end
    spBm = rxnBiomassID;
end
%ATPM
if iscell(spATPM)
    spATPM0 = spATPM;
    spATPM = zeros(nSp,1);
    for j = 1: nSp
        spATPM(j) = findRxnIDs(modelCell{j}, spATPM0{j});
    end
end
if isempty(addExRxns)
    addExRxns = true;
end
if isempty(sepUtEx)
    sepUtEx = true;
end
msg = '';
%% Copy fields from COBRA model
field = {};
[fieldNumeric, fieldCell, fieldStruct]  = deal(false(0));
for jSp = 1:nSp
    fCur = fieldnames(modelCell{jSp});
    fNcur = cellfun(@(x) isnumeric(modelCell{jSp}.(x)),fCur);
    fCcur = cellfun(@(x) iscell(modelCell{jSp}.(x)),fCur);
    fScur = cellfun(@(x) isstruct(modelCell{jSp}.(x)),fCur);
    field = [field; fCur];
    fieldNumeric = [fieldNumeric; fNcur];
    fieldCell = [fieldCell; fCcur];
    fieldStruct = [fieldStruct; fScur];
end
[field,id] = unique(field);
[fieldNumeric,fieldCell,fieldStruct] = deal(fieldNumeric(id),fieldCell(id),fieldStruct(id));
%fields that need special care
id = ismember(field, {'S', 'rxns', 'mets', 'rev', 'lb', 'ub', 'c', 'b', ...
    'metComs','rxnGeneMat','genes', 'grRules', 'rules', 'comps', 'compNames'});
[field,fieldNumeric,fieldCell,fieldStruct] = deal(field(~id),fieldNumeric(~id),fieldCell(~id),fieldStruct(~id));
rxnField = unique([rxnField(:);{'confidenceScores';'subSystems'}]);

modelCom = struct();
for j = 1:numel(field)
    modelCom.(field{j}) = [];
end

%% fields to be changed
[S, rxns, mets, metsCom, rev, lb, ub, c, b] = deal([]);
%map rxns, mets, genes, compartments to species
[rxnSps, metSps, geneSps, compSps] = deal([]);
%rxn IDs of biomass
spBmId = zeros(nSp, 1);

%extracellular metabolites
ex = cell(nSp, 1);
%exchange reactions
rxnEx = cell(nSp, 1);
%map ex mets and ex rxns
rxnEx2met = cell(nSp, 1);

% compartments
[comps, compNames] = deal({});

%% loop for each species
if ~isempty(metField)
    metFieldKnown = ismember(field,metField);
else
    metFieldKnown = false(numel(field),1);
end
if ~isempty(rxnField)
    rxnFieldKnown = ismember(field,rxnField);
else
    rxnFieldKnown = false(numel(field),1);
end
if ~isempty(geneField)
    geneFieldKnown = ismember(field, geneField);
else
    geneFieldKnown = false(numel(field),1);
end
if ~isempty(compField)
    compFieldKnown = ismember(field, compField);
else
    compFieldKnown = false(numel(field),1);
end
%met-related fields
metFieldL = strncmp('met',field,3) | metFieldKnown | ismember(field, {'csense'});
%rxn-related fields
rxnFieldL = strncmp('rxn',field,3) | rxnFieldKnown;
% gene-related fields
geneFieldL = strncmp('gene',field,4) | ismember(field, {'proteins'}) | geneFieldKnown;
% compartment-related fields
compFieldL = strncmp('comp',field,4) | compFieldKnown;
%all other fields just put in a cell for each model
spFieldL = ~metFieldL & ~rxnFieldL & ~geneFieldL & ~compFieldL;

msg = printAndRecordMsg(msg, sprintf('Compile metabolites and reactions...\n'), verbFlag);

row = 0;
col = 0;
for j = 1:nSp
    [rowJ, colJ] = size(modelCell{j}.S);
    if ~isfield(modelCell{j},'genes')
        modelCell{j}.genes = {};
    end
    geneJ = numel(modelCell{j}.genes);
    modelMsg = printCurrentModel(j, spAbbr, spName);
    if isfield(modelCell{j}, 'metComs') %given the mapping to community metabolites
        if numel(modelCell{j}.metComs) < numel(modelCell{j}.mets)
            msg = warnAndRecordMsg(msg, sprintf('%s: size of metComs < size of mets.', modelMsg), verbFlag);
            modelCell{j}.metComs(end+1:numel(modelCell{j}.mets)) = {''};
        elseif numel(modelCell{j}.metComs) > numel(modelCell{j}.mets)
            msg = warnAndRecordMsg(msg, sprintf('%s: size of metComs > size of mets.', modelMsg), verbFlag);
            modelCell{j}.metComs(numel(modelCell{j}.mets)+1:end) = [];
        end
        %logical vector for extracellular metabolites
        ex{j} = ~cellfun(@isempty, modelCell{j}.metComs);
        modelCell{j}.metComs(ex{j}) = regexprep(modelCell{j}.metComs(ex{j}),'\[[^\[\]]*\]$','');
    else
        if ~isempty(metExId) %using identifier for extracellular mets
            %logical vector for extracellular metabolites
            ex{j} = cellfun(@(x) ~isempty(strfind(x, metExId)), modelCell{j}.mets);
        else %if no identifier, just check exchange reactions (will include sink/demand reactions)
            ex{j} = any(modelCell{j}.S(:, sum(modelCell{j}.S ~= 0) == 1), 2);
        end
        modelCell{j}.metComs = repmat({''},rowJ,1);
        modelCell{j}.metComs(ex{j}) = regexprep(modelCell{j}.mets(ex{j}),'\[[^\[\]]*\]$','');
    end
    %logical vector for exchange reactions
    rxnEx{j} = (sum(modelCell{j}.S(ex{j},:) ~= 0) == 1)' & (sum(modelCell{j}.S(~ex{j},:) ~= 0) == 0)';
    if addExRxns
        % add exchange reactions for those mets in ex{j} having no exchange rxns in rxnEx{j}
        metExWoER = find(ex{j} & ~any(modelCell{j}.S(:, rxnEx{j}), 2));
        if ~isempty(metExWoER)
            msg = printAndRecordMsg(msg, sprintf('%s:\n', modelMsg), verbFlag);
            for mJ = 1:numel(metExWoER)
                rxn = ['EX_' regexprep(modelCell{j}.mets{metExWoER(mJ)}, '\[(.*)\]', '\($1\)')];
                modelCell{j} = addReaction(modelCell{j}, rxn, ...
                    'reactionName', ['Exchange reaction for ' modelCell{j}.metNames{metExWoER(mJ)}], ...
                    'reactionFormula', [modelCell{j}.mets{metExWoER(mJ)}, ' <=>'], ...
                    'printLevel', 0);
                msg = printAndRecordMsg(msg, sprintf('\tExchang reaction added for %s, %s: %s\n', ...
                    modelCell{j}.mets{metExWoER(mJ)}, rxn, [modelCell{j}.mets{metExWoER(mJ)}, ' <=>']), verbFlag);
            end
            % redefine the number of reactions and exchange reaction indices
            colJ = colJ + numel(metExWoER);
            rxnEx{j} = (sum(modelCell{j}.S(ex{j},:) ~= 0) == 1)' & (sum(modelCell{j}.S(~ex{j},:) ~= 0) == 0)';
            % update rxn fields
            for jF = 1:numel(field)
                if rxnFieldL(jF) && size(modelCell{j}.(field{jF}), 1) < colJ
                    if isnumeric(modelCell{j}.(field{jF}))
                        modelCell{j}.(field{jF})((end + 1):colJ, :) = 0;
                    else
                        modelCell{j}.(field{jF})((end + 1):colJ, :) = {''};
                    end
                end
            end
        end
    end
    %get bounds and objective
    lbJ = modelCell{j}.lb;
    ubJ = modelCell{j}.ub;
    cJ = modelCell{j}.c;
    %indices for exchange reactions mapped to metabolites
    rxnEx2met{j} = zeros(rowJ,2);
    for k = 1:colJ
        if rxnEx{j}(k)
            metJK = find(modelCell{j}.S(:, k), 1);
            rxnEx2met{j}(metJK,:) = [col + k, k];
            %set positive flux of exchange reaction as uptake
            %(convention)
            if modelCell{j}.S(metJK, k) > 0
                s = modelCell{j}.S(metJK,k);
                modelCell{j}.S(metJK,k) = -1;
                [ubJ(k), lbJ(k), cJ(k)] = deal(-s * lbJ(k), -s * ubJ(k), cJ(k));
            end
        end
    end
    %update community metabolites (now the model has .metComs, which has no compartment identifier)
    metsCom = unique([metsCom; modelCell{j}.metComs(ex{j})]);    
    %stoichiometric matrix    
    S = [S                 sparse(row, colJ);...
         sparse(rowJ, col) sparse(modelCell{j}.S)];
    %reversibility, bounds, objective and RHS
    if ~isfield(modelCell{j}, 'rev')
        modelCell{j}.rev = double(modelCell{j}.lb < 0);
    end
    rev = [rev; modelCell{j}.rev];
    lb = [lb; lbJ];
    ub = [ub; ubJ];
    c = [c; cJ];
    b = [b; modelCell{j}.b];
    %add species's name to rxns and mets (add to the compartment if exist)
    metJ = regexprep(modelCell{j}.mets,'\]$',['_' spAbbr{j} '\]']);
    metWoComp = cellfun(@(x) ~strcmp(x(end),']'), metJ);
    metJ(metWoComp) = strcat(metJ(metWoComp), '[', spAbbr{j}, ']');
    mets = [mets; metJ];
    rxns = [rxns; strcat(modelCell{j}.rxns,'_',spAbbr{j})];
    
    % compartments
    if isfield(modelCell{j}, 'comps')
        % use the compartment IDs
        compsJ = modelCell{j}.comps;
    else
        % retreive the compartment IDs from .mets if not exist
        compsJ = regexp(modelCell{j}.mets, '(\[.*\])$','tokens','once');
        % if there is any met without compartment ID
        emptyComp = cellfun(@isempty, compsJ);
        if any(emptyComp)
            compsJ = compsJ(~emptyComp);
            emptyComp = true;
        else
            emptyComp = false;
        end
        compsJ = cellfun(@(x) x{1}, compsJ ,'UniformOutput',false);
        compsJ = unique(compsJ(:));
        if emptyComp
            compsJ{end + 1} = '';
        end
    end
    comps = [comps; strcat(compsJ, '_', spAbbr{j})];
    if isfield(modelCell{j}, 'compNames')
        compNames  = [compNames; strcat(modelCell{j}.compNames, '_', spName{j})];
    else
        compNames  = [compNames; strcat(compsJ, '_', spName{j})];
    end
    nCompJ = numel(compsJ);
    
    %incorporate other field
    for k = 1:numel(field)
        if metFieldL(k)
            str = 'mets';
            sizeCk = rowJ;
        elseif rxnFieldL(k)
            str = 'rxns';
            sizeCk = colJ;
        elseif geneFieldL(k)
            str = 'genes';
            sizeCk = geneJ;
        elseif compFieldL(k)
            str = 'comps';
            sizeCk = nCompJ;
        else
            %spField
            str = '1';
            sizeCk = 1;
        end
        if isfield(modelCell{j}, field{k})
            if spFieldL(k)
                % Whatever it is, put it in a cell.
                fieldJK = {modelCell{j}.(field{k})};
                modelCom.(field{k}) = [modelCom.(field{k}); fieldJK];
            else
                %check sizes
                if size(modelCell{j}.(field{k}),1) ~= sizeCk
                    if size(modelCell{j}.(field{k}),2) == sizeCk
                        fieldJK = modelCell{j}.(field{k})';
                    else
                        error('Dimension of modelCell{%d}.%s (%d,%d) does not match %s (%d).',...
                            j,field{k},size(modelCell{j}.(field{k})),str,sizeCk)
                    end
                else
                    fieldJK = modelCell{j}.(field{k});
                end
                %assignment
                if isempty(modelCom.(field{k}))
                    modelCom.(field{k}) = fieldJK;
                else
                    modelCom.(field{k}) = [modelCom.(field{k});fieldJK];
                end
            end
        else
            if spFieldL(k) && isempty(modelCom.(field{k}))
                modelCom.(field{k}) = {[]};
            elseif spFieldL(k) && ~isempty(modelCom.(field{k}))
                modelCom.(field{k})(end + 1, 1) = {[]};
            else
                if isempty(modelCom.(field{k}))
                    if fieldNumeric(k)
                        modelCom.(field{k}) = zeros(sizeCk,1);
                    elseif fieldCell(k)
                        modelCom.(field{k}) = repmat({''},sizeCk,1);
                    elseif fieldStruct(k)
                        modelCom.(field{k}) = repmat(struct(),sizeCk,1);
                    end
                else
                    if fieldNumeric(k)
                        modelCom.(field{k})(end+1:end+sizeCk,:) = 0;
                    elseif fieldCell(k)
                        modelCom.(field{k})(end+1:end+sizeCk,:) = {''};
                    elseif fieldStruct(k)
                        f = fieldnames(modelCom.(field{k})(end));
                        modelCom.(field{k})(end+1:end+sizeCk,:) = ...
                            repmat(cell2struct(repmat({[]},numel(f),1),f),sizeCk,1);
                    end
                end
            end
        end
    end
    
    %species specific rxns and mets
    rxnSps = [rxnSps; j * ones(colJ,1)];
    metSps = [metSps; j * ones(rowJ,1)];
    geneSps = [geneSps; j * ones(geneJ, 1)];
    compSps = [compSps; j * ones(nCompJ, 1)];
    %biomass rxn ID
    spBmId(j) = col + spBm(j);
    if ~isempty(spATPM)
        spATPM(j) = col + spATPM(j);
    end
    %size of the network
    row = row + rowJ;
    col = col + colJ;
end

%% Community metabolites
msg = printAndRecordMsg(msg, sprintf('Add community metabolites and exchange reactions...\n'), verbFlag);
metsCom = sort(unique(metsCom));
%Ids of exchange reactions corresponding to community metabolites
% [a_ij] = exchange reaction Id for community metabolite i and species j 
EXsp = zeros(numel(metsCom), nSp);
Msp = zeros(numel(metsCom), nSp);
[rowS,colS,entryS] = deal([]);
for kSp = 1:nSp
    %organism-community exchange reactions
    [r0,c0,e0] = find(modelCell{kSp}.S);
    modelCell{kSp}.metComs(cellfun(@isempty,modelCell{kSp}.metComs)) = {''};
    [yn,id] = ismember(modelCell{kSp}.metComs,metsCom);
    % in case some extracellular metabolites do not have exchange
    % reactions. Then organism-community exchange reactions are not added.
    ynCom = yn & rxnEx2met{kSp}(:, 2) > 0;
    rowS = [rowS; id(ynCom)];
    colS = [colS; rxnEx2met{kSp}(ynCom,1)];
    [yn2,id2] = ismember([find(ynCom),rxnEx2met{kSp}(ynCom,2)],[r0,c0],'rows');
    entryS = [entryS; -e0(id2)];
    EXsp(id(yn),kSp) = rxnEx2met{kSp}(yn,1);
    Msp(id(yn),kSp) = find(metSps == kSp, 1) - 1 + find(yn);
    metNotAdded = find(yn & ~ynCom);
    if ~isempty(metNotAdded)
        modelMsg = printCurrentModel(kSp, spAbbr, spName);
        msg = printAndRecordMsg(msg, sprintf(['%s\n\tThe following mets are not added as community mets', ...
            'because there are no exchange reactions for them in the model:\n'], modelMsg), verbFlag);
        for mK = 1:numel(metNotAdded)
            msg = printAndRecordMsg(msg, sprintf('\t%s\t%s\n', modelCell{kSp}.mets{metNotAdded(mK)}, ...
                modelCell{kSp}.metComs{metNotAdded(mK)}), verbFlag);
        end
    end
end
% community exchange reactions are NOT added for extracellular metabolites
% without any exchange reactions in any organisms
wtExCom = any(EXsp, 2);
mCom = sum(wtExCom);
[~, rowS] = ismember(rowS, find(wtExCom)); 
if sepUtEx
    %community uptake/export reactions
    rowS = [rowS; repmat((1:mCom)',2,1)];
    colS = [colS; ((col + 1):(col + 2 * mCom))'];
    entryS = [entryS; ones(mCom, 1); -ones(mCom, 1)];
    SmetCom = sparse(rowS, colS, entryS, mCom, col + 2 * mCom);
    S = [S sparse([],[],[], row, 2 * mCom); SmetCom];
    % default no uptake
    %For non-limiting substrate for uptake, if the lower bound for uptake is high (say
    % 1000), then the upper bound of the corresponding export reaction of the
    % community metabolites should be significantly larger (say 10000) in order
    % not to overconstrain the community in the way that it is not allowed to
    % produce those metabolites
    ub = [ub; zeros(mCom, 1); 10000 * ones(mCom, 1)];
    lb = [lb; zeros(2 * mCom, 1)];
    c = [c;  zeros(2 * mCom, 1)];
    rev = [rev;  zeros(2 * mCom, 1)];
    rxnSps = [rxnSps; zeros(2 * mCom, 1)];
    %names of community exchange reactions
    rxns = [rxns; strcat('UT_',metsCom(wtExCom),'(u)'); strcat('EX_',metsCom(wtExCom),'(u)')];
else
    rowS = [rowS; (1:mCom)'];
    colS = [colS; ((col + 1):(col + mCom))'];
    entryS = [entryS; -ones(mCom, 1)];
    SmetCom = sparse(rowS, colS, entryS, mCom, col + mCom);
    S = [S sparse([],[],[], row, mCom); SmetCom];
    ub = [ub; 1000 * ones(mCom, 1)];
    % default no uptake
    lb = [lb; zeros(mCom, 1)];
    c = [c;  zeros(mCom, 1)];
    rev = [rev;  ones(mCom, 1)];
    rxnSps = [rxnSps; zeros(mCom, 1)];
    %names of community exchange reactions
    rxns = [rxns; strcat('EX_',metsCom(wtExCom),'(u)')];
end
b = [b; zeros(mCom, 1)];
metSps = [metSps; zeros(mCom, 1)];
%names of community metabolites
mets = [mets; strcat(metsCom(wtExCom),'[u]')];

comps = [comps; {'u'}];
compNames = [compNames; {'Community'}];
compSps(end + 1) = 0;

[modelCom.rxns, modelCom.mets, modelCom.S, modelCom.c, modelCom.lb, ...
    modelCom.ub, modelCom.b, modelCom.rev, modelCom.comps, modelCom.compNames] =...
    deal(rxns, mets, S, c, lb, ub, b, rev, comps, compNames);
%get community reaction indices
indCom = struct();

rxnSD = sum(modelCom.S ~= 0, 1) <= 1;
rxnSD((col + 1) : (col + end)) = false;
rxnSD = find(rxnSD);
%reaction Ids for [uptake | export] of community metabolites
if sepUtEx
    EXcom = [((col + 1): (col + mCom))', ((col + mCom + 1): (col + 2 * mCom))'];
else
    EXcom = (col + 1: col + mCom)';
end
EXsp = EXsp(wtExCom, :);
[indCom.spBm, indCom.spATPM, indCom.rxnSD, indCom.EXcom, indCom.EXsp,...
    indCom.Mcom, indCom.Msp, indCom.rxnSps, indCom.metSps, indCom.geneSps, indCom.compSps] ...
    = deal(spBmId, spATPM, rxnSD, EXcom, EXsp, ...
    ((row + 1) : (row + mCom))', Msp(wtExCom, :), rxnSps, metSps, geneSps, compSps);


%add rxnNames if exist                
if ~isfield(modelCom, 'rxnNames')
    modelCom.rxnNames = modelCom.rxns;
else
    if numel(modelCom.rxnNames) ~= col
        if sepUtEx
            modelCom.rxnNames((col + 1): (col + 2 * mCom)) = modelCom.rxns((col + 1): (col + 2 * mCom));
        else
            modelCom.rxnNames((col + 1): (col + mCom)) = modelCom.rxns((col + 1): (col + mCom));
        end
    else
        for j = 1:mCom
            rxnNameLength = zeros(nSp, 1);
            for k = 1:nSp
                if EXsp(j,k) > 0
                    rxnNameLength(k) = length(modelCom.rxnNames{EXsp(j,k)});
                end
            end
            [maxLength, maxLengthId] = max(rxnNameLength);
            if sepUtEx
                if maxLength > 0
                    rxnNameJ = modelCom.rxnNames{EXsp(j,maxLengthId)};
                    if ~isempty(strfind(rxnNameJ, 'exchange'))
                        rxnNameJut = strrep(rxnNameJ, 'exchange', '(community uptake)');
                        rxnNameJex = strrep(rxnNameJ, 'exchange', '(community export)');
                    else
                        rxnNameJut = strcat(rxnNameJ, ' (community uptake)');
                        rxnNameJex = strcat(rxnNameJ, ' (community export)');
                    end
                    modelCom.rxnNames{col + j} = rxnNameJut;
                    modelCom.rxnNames{col + mCom + j} = rxnNameJex;
                else
                    modelCom.rxnNames{col + j} = modelCom.rxns{col + j};
                    modelCom.rxnNames{col + mCom + j} = modelCom.rxns{col + mCom + j};
                end
            else
                if maxLength > 0
                    rxnNameJ = modelCom.rxnNames{EXsp(j,maxLengthId)};
                    if ~isempty(strfind(rxnNameJ, 'exchange'))
                        rxnNameJex = strrep(rxnNameJ, 'exchange', '(community export)');
                    else
                        rxnNameJex = strcat(rxnNameJ, ' (community export)');
                    end
                    modelCom.rxnNames{col + j} = rxnNameJex;
                else
                    modelCom.rxnNames{col + j} = modelCom.rxns{col + j};
                end
            end
        end
    end
end
%add metNames, metFormulas and other metField if exist
if ~isfield(modelCom, 'metFormulas')
    modelCom.metFormulas = repmat({''}, row + mCom, 1);
end
if ~isfield(modelCom, 'metNames')
    modelCom.metNames = modelCom.mets;
else
    if numel(modelCom.metNames) ~= row
        modelCom.metNames((row + 1): (row + mCom)) = modelCom.mets((row + 1): (row + mCom));
    else
        for j = 1:mCom
%             metNameLength = zeros(nSp, 1);
            metNameId = zeros(nSp, 1);
            metNameJ = {};
            metForm = '';
            %get all names, put in cell array
            for k = 1:nSp
                if EXsp(j,k) > 0
                    metEUid = modelCom.S(:,EXsp(j,k)) ~= 0;
                    metEUid(row + j) = false;
                    metNameId(k) = find(metEUid, 1);
                    metNameJK = modelCom.metNames{metNameId(k)};
                    metNameJK = strrep(strrep(metNameJK, '(extracellular)', ''), '(Extracellular)', '');
                    metNameJK = strrep(strrep(metNameJK, 'extracellular', ''), 'Extracellular', '');
                    metNameJK = unique(strtrim(splitString(strtrim(metNameJK),';')))';
                    metNameJ = [metNameJ metNameJK];
                    if isfield(modelCom, 'metFormulas') && isempty(metForm)
                        metForm = strtrim(modelCom.metFormulas{metNameId(k)});
                    end
                end
            end
            modelCom.metFormulas{row + j} = metForm;
            metNameJ = unique(metNameJ);
            %remove those are contained totally in another name
            %             k = 1;
            %             while k <= numel(metNameJ)
            %                 if any(cellfun(@(x) ~isempty(strfind(x, lower(metNameJ{k}))), lower(metNameJ([1:k-1 k+1:end]))))
            %                     %if totally contained, delete it
            %                     metNameJ(k) = [];
            %                 else
            %                     k = k + 1;
            %                 end
            %             end
            metNameLength = cellfun(@length,metNameJ);
            if ~isempty(metNameLength)
                metNameJ = metNameJ(metNameLength > 0);
            end
            if ~isempty(metNameJ)
                if iscell(metNameJ) && numel(metNameJ) == 1
                    modelCom.metNames{row + j} = metNameJ{1};
                else
                    modelCom.metNames{row + j} = strjoin(metNameJ,'|');
                end
            else
                modelCom.metNames{row + j} = modelCom.mets{row + j};
            end
            for kF = 1:numel(metField)
                modelCom.(metField{kF})(row + j) = modelCom.(metField{kF})(metEUid);
            end
                
%             for k = 1:nSp
%                 if modelCom.EXsp(j,k) > 0
%                     metEUid = modelCom.S(:,modelCom.EXsp(j,k)) ~= 0;
%                     metEUid(row + j) = false;
%                     metNameId(k) = find(metEUid, 1);
%                     metNameLength(k) = length(modelCom.metNames{metNameId(k)});
%                 end
%             end
%             [maxLength, maxLengthId] = max(metNameLength);
%             if maxLength > 0
%                 modelCom.metNames{row + j} = modelCom.metNames{metNameId(maxLengthId)};
%             else
%                 modelCom.metNames{row + j} = modelCom.mets{row + j};
%             end
        end
    end
end
if isfield(modelCom,'subSystems')
    if sepUtEx
        modelCom.subSystems((col + 1): (col + mCom * 2)) = repmat({'community exchange'}, 2 * mCom, 1);
    else
        modelCom.subSystems((col + 1): (col + mCom)) = repmat({'community exchange'}, mCom, 1);
    end
end
%extend all other fields
uMet = zeros(mCom,1);
Msp = Msp(wtExCom, :);
for jM = 1:mCom
    uMet(jM) = Msp(jM, find(Msp(jM, :), 1));
end
for j = 1:numel(field)
    if metFieldL(j) || rxnFieldL(j)
        if metFieldL(j) && ~any(strcmp({'metNames', 'metFormulas'}, field{j}))
            v = (row + 1): (row + mCom);
            modelCom.(field{j})(v, :) = modelCom.(field{j})(uMet, :);
        elseif rxnFieldL(j) && ~any(strcmp({'subSystems', 'rxnNames'}, field{j}))
            if sepUtEx
                v = (col + 1): (col + mCom * 2);
            else
                v = (col + 1): (col + mCom);
            end
            if fieldNumeric(j)
                u = 0;
            elseif fieldCell(j)
                u = {''};
            elseif fieldStruct(j)
                f = fieldnames(modelCom.(field{j}));
                u = cell2struct(repmat({[]},numel(f),1),f);
            end
            modelCom.(field{j})(v, :) = u;
        end
    end
end

%% special care for genes, rxnGeneMat, grRules and rules
msg = printAndRecordMsg(msg, sprintf('Merge gene-reaction rules...\n'), verbFlag);
[rGMrow,rGMcol,rGMent] = deal([]);
rGMm = 0;
rGMn = 0;
[genes, grRules] = deal({});
[modelCom.grRules, modelCom.rules] = deal(repmat({''}, numel(modelCom.rxns), 1));
modelCom.genes = {};
for jG = 1:numel(geneField)
    modelCom.(geneField{jG}) = {};
end
for jSp = 1:nSp
    modelMsg = printCurrentModel(jSp, spAbbr, spName);
    if ~isfield(modelCell{jSp}, 'rules') 
        if isfield(modelCell{jSp}, 'grRules')
            modelCell{jSp} = generateRules(model, 0);
        else
            modelCell{jSp}.rules = repmat({''}, numel(modelCell{jSp}.rxns), 1);
        end
    end
    % add organism's tag to genes to ensure unique gene names in the model
    modelCell{jSp}.genes = strcat(modelCell{jSp}.genes(:), '_', spAbbr{jSp});
    % update grRules
    modelCell{jSp} = creategrRulesField(modelCell{jSp});
    
    genes = [genes; modelCell{jSp}.genes];
    grRules = [grRules; modelCell{jSp}.grRules];
end
modelCom.genes = genes;
% no grRules for community exchange reactions
if sepUtEx
    grRules = [grRules; repmat({''}, mCom * 2, 1)];
else
    grRules = [grRules; repmat({''}, mCom, 1)];
end
modelCom.grRules = grRules;
% build .rules from .grRules
modelCom = generateRules(modelCom, 0);
% build .rxnGeneMat
modelCom = buildRxnGeneMat(modelCom);

% get community info
infoCom = infoCom2indCom(modelCom,indCom,true,spAbbr,spName);
% add infoCom and indCom into modelCom
modelCom.infoCom = infoCom;
modelCom.indCom = indCom;

if isfield(modelCom, 'description') && iscell(modelCom.description)
    modelCom.description = ['Community model built from ' strjoin(modelCom.description, ', ')];
end
if isfield(modelCom, 'osenseStr')
    modelCom.osenseStr = 'max';
end
end

function str = printCurrentModel(kSp, spAbbr, spName)
str = sprintf('Input model %d %s', kSp, spAbbr{kSp});
if ~strcmp(spName{kSp}, spAbbr{kSp})
    str = sprintf('%s\t%s', str, spName{kSp});
end
end

function msg = printAndRecordMsg(msg, msgCur, verbFlag)
if verbFlag, fprintf(msgCur); end
msg = sprintf('%s%s', msg, msgCur);
end

function msg = warnAndRecordMsg(msg, msgCur, verbFlag)
if verbFlag, warning(msgCur); end
msg = sprintf('%s%s', msg, msgCur);
end

