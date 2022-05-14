# minMicrobiome
minMicrobiome is an algorithm to identify minimal microbiomes that are capable of a given functionality (metabolite production) from a large community of microbes, using a constraint-based approach. To arrive at a minimal microbiome, sequential deletion of member species is done and is followed by solving a mixed integer linear programming problem (MILP) with the objective to minimize L1 norm of the membership vector. minMicrobiome function accepts the type of functionality to be retained as one of the input parameters. The identified minimal microbiome need not be unique and running multiple iterations of the algorithm can yield more possible minimal microbiomes.

## Requirements
1. COBRA toolbox
2. Solver for MILP and LP such as Gurobi

## Instructions for use
The function can be called using the command:
`[solutionPert, supp_out] = MinMicrobiome(modelCom, options);`

### Inputs
_modelCom_ - a community model of AGORA models created using the COBRA toolbox function `createCommModel`, such that:  
- Each species in the community is represented as ‘_ _org1_’, ‘_ _org2_’, etc.   
- Exchange reactions start with ‘_EX_ _’  
- Names of biomass equations start with ‘_biomass_’  

### Optional inputs
1.	Fraction of growth rate to be retained in minimal microbiome – _gr_frac_ (default value = 0.8)  
2.	Fraction of SCFA to be retained in minimal microbiome – _scfa_frac_ (default value = 0.8)
3.	For coupling biomass reaction to the others in a species: coupling constraint, _c_ (default value = 1000) and threshold _u_ (default value = 0.01)
4.	The functionality of the community that needs to be preserved – _constraint_, namely:
  - 1: maximum weighted sum of SCFAs (metabolites 1, 2 and 3, default _wt_ = 1:1:1),
  - 2: maximum metabolite 1 production 
  - 3: maximum metabolite 2 production 
  - 4: maximum metabolite 3 production     
  (default _constraint_ value = 1. Default _met_names_: metabolite 1 = acetate, metabolite 2 = butyrate and metabolite 3 = propionate)      
 5.	Maximum number of organisms for MILP – _maxMILP_ (default value = 10)  
 6.	Total no. of iterations the code should run – _iter_ (default value = 3)  
 7.	Maximum number of rounds the generated random sequence should be ran over for deletion (MILP may not reduce to the desired number with 1 round) – _max_del_rounds_ (default value = 10)   
 8.	Fraction of original growth rates for constraint while finding maximum SCFA production of the original community - _gr_opt_frac_ (default value = 0.99)  
 9. Sequence for deletion - _del_seq_ (generated randomly by default)  

### Outputs
1.  _solutionPert_: a cell containing the identified minimal microbiomes after running the specified number of iterations (No. of minimal microbiomes calculated <= no. of iterations) 
2.  _supp_out_: the additional output containing the full community growth rate (_gr_max_), maximum individual growth rates and maximum SCFA production (_max_scfa_) of the original community, and the list of organisms removed during the initial deletion 

