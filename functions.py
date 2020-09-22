

def join_lumen_rxns(model,mapping,Metabolite,Reaction):
    #Makes lumen exchange reactions and add them to the model
    
    for i in range(len(model.reactions)):
        
        if model.reactions[i].id.startswith('Ex') and len(model.reactions[i].metabolites)==1:

            model.reactions[i].lower_bound = 0
            
            if model.reactions[i].id[3:] in mapping['coreco'].values:
                
                s = mapping[mapping['coreco']==model.reactions[i].id[3:]].iloc[0,0]
                new_met = Metabolite(s+'_u',compartment = 'u')

                #add the lumen metabolite to the right side of the exchange reaction
                model.reactions[i].add_metabolites({new_met:1})

                model.reactions[i].bounds = (-1000,1000)

                #add an exchange reaction for the lumen metabolite (lumen <=> outside)
                r = Reaction("EX_u_"+s)
                model.add_reaction(r)
                r.add_metabolites({new_met:-1})
            
    return model
            


def constrain_model(model,diet,identifier):
    #Constraints the model based on the given diet
    for i in range(len(model.reactions)):
        if model.reactions[i].id[:len(identifier)]==identifier:
            
            if model.reactions[i].id in diet.keys():
                model.reactions[i].lower_bound = diet[model.reactions[i].id]
            else:
                model.reactions[i].lower_bound=0
            
    return model



   
def change_ids(model):
    #Changes metabolite and reaction IDs of the bacterial model so that it is distinguishable.     
    for i in range(len(model.reactions)):
        model.reactions[i].id = 'model2_' + model.reactions[i].id

    for i in range(len(model.metabolites)):
        model.metabolites[i].id = 'model2_' + model.metabolites[i].id

    return model



def tune_bac_lumen_rxns(bac_model,calb_model,Metabolite,Reaction,source):
    #Creates corresponding lumen reactions to prepare the bacterial model to be joined with Candida albicans model
    second_lumen_ext_rxns = list()
    second_lumen_ext_mets = list()
    
    for i in range(len(bac_model.reactions)):
        rxn_id = bac_model.reactions[i].id
        
        if rxn_id.startswith('model2_EX_'):
            
            bac_model.reactions[i].bounds = (-1000,1000)
            
            lumen_met = list(bac_model.reactions[i].metabolites)[0].id
            lumen_met = lumen_met.replace('model2_','')
            
            if source == 'agora':
                lumen_met = lumen_met.replace('[e]','_u')
            elif source == 'carveme':
                lumen_met = lumen_met.replace('_e','_u')

            #check if the id comes from carveme model that may contain "__" in the id. In this case, '__' should be replaced with '_' in the id
            if '__' in lumen_met:
                lumen_met = lumen_met.replace('__','_')

            new_met = Metabolite(lumen_met , compartment = 'u')
            bac_model.reactions[i].add_metabolites({new_met:1})
            
            if lumen_met not in calb_model.metabolites:
                
                if source == 'agora':
                    new_rxn_id = 'EX_u_' + bac_model.reactions[i].id[10:-3]
                    
                elif source == 'carveme':
                    new_rxn_id = 'EX_u_' + bac_model.reactions[i].id[10:-2]
                    if '__' in new_rxn_id:
                        new_rxn_id = new_rxn_id.replace('__','_')
                        
                new_rxn = Reaction(new_rxn_id)

                second_lumen_ext_rxns.append(new_rxn)
                second_lumen_ext_mets.append(new_met)

    return bac_model,second_lumen_ext_rxns,second_lumen_ext_mets
    

def make_joint_model(calb_model,bac_model,second_lumen_ext_rxns,second_lumen_ext_mets,diet,copy,source):
    #Creates the joint model
    joint_model = copy.deepcopy(calb_model)
    joint_model.solver = 'cplex'
    
    
    for i in range(len(second_lumen_ext_rxns)):
            
        joint_model.add_reactions([second_lumen_ext_rxns[i]])
            
        second_lumen_ext_rxns[i].add_metabolites({second_lumen_ext_mets[i]:-1})
        second_lumen_ext_rxns[i].bounds = (0,1000)
            
        if (second_lumen_ext_rxns[i].id in diet.keys()):
            second_lumen_ext_rxns[i].lower_bound = diet[second_lumen_ext_rxns[i].id]


    joint_model.add_reactions(bac_model.reactions)

    #Set for anaerobic growth
    o2 = joint_model.reactions.get_by_id('EX_u_o2')
    o2.lower_bound = 0

    if source == 'agora':
        bac_biomass_id = 'model2_biomass'
    elif source == 'carveme':
        bac_biomass_id = 'model2_Growth'
    
    for i in range(len(joint_model.reactions)):
        if joint_model.reactions[i].id.startswith(bac_biomass_id):
            biomass_id_second_model = joint_model.reactions[i].id    
            
    biomass1 = joint_model.reactions.get_by_id('scerbiomasspseudoreaction')
    biomass2 = joint_model.reactions.get_by_id(biomass_id_second_model)

    biomass1.lower_bound = 0.001
    biomass2.lower_bound = 0.001
    
    joint_model.objective = [biomass1,biomass2]

    return joint_model


def apply_coupling_constraints(cpx,cplex,source):
    #Adds the coupling constrains to the joint model
    vars_num = cpx.variables.get_num()

    if source == 'agora':
        bac_biomass_id = 'model2_biomass'
    elif source == 'carveme':
        bac_biomass_id = 'model2_Growth'
    
    #make a loop to find the location of biomass reactions for both models in cpx
    for i in range(vars_num):
        
        ss = cpx.variables.get_names(i)
        
        if ss.startswith('scerbiomasspseudoreaction') and ss.find('reverse')==-1:
            s = cpx.variables.get_names(i)
            loc_bio_1 = i
        
        if ss.startswith(bac_biomass_id) and ss.find('reverse')==-1:
            loc_bio_2 = i
        


    #take all the variables(reactions) except exchanges and biomass production
    first_forward_vars = list()
    first_backward_vars = list()
    second_forward_vars = list()
    second_backward_vars = list()

    vars_num = cpx.variables.get_num()

    bac_biomass_id = bac_biomass_id.replace('model2_','')
    
    for i in range(vars_num):
        
        s = cpx.variables.get_names(i)
        
        if not s.startswith('_EX_u_') and \
        not s.startswith('_Ex_') and \
        not s.startswith('model2_EX_') and \
        bac_biomass_id not in s and \
        not s.startswith('scerbiomasspseudoreaction'):
            
            if not s.startswith('model2'):

                if i%2 == 0:
                    first_forward_vars.append(s)
                else:
                    first_backward_vars.append(s)

            else:
                if i%2 == 0:
                    second_forward_vars.append(s)
                else:
                    second_backward_vars.append(s)
        


    #add biomass reaction at the end of the variables' list
    first_forward_vars.append(cpx.variables.get_names(loc_bio_1))
    first_backward_vars.append(cpx.variables.get_names(loc_bio_1))
    second_forward_vars.append(cpx.variables.get_names(loc_bio_2))
    second_backward_vars.append(cpx.variables.get_names(loc_bio_2))
    
    #couple the candida albicans reactions to the corresponding biomass reaction 
    i=0
    while i<len(first_forward_vars)-1:
        
        cpx.linear_constraints.add(
            lin_expr = [cplex.SparsePair(ind = [first_forward_vars[i],first_forward_vars[len(first_forward_vars)-1]],val = [1, -400]),
                        cplex.SparsePair(ind = [first_backward_vars[i],first_backward_vars[len(first_forward_vars)-1]],val = [1, -400])],
            
            senses = ["L","L"],
            rhs=[0.01,0.01],
            range_values = [0,0],
            names=["coupled_"+first_forward_vars[i], "coupled_"+first_backward_vars[i]])

        i += 1


    #couple the bacteria reactions to the corresponding biomass reaction 
    i=0
    while i<len(second_forward_vars)-1:
        cpx.linear_constraints.add(
            lin_expr = [cplex.SparsePair(ind = [second_forward_vars[i],second_forward_vars[len(second_forward_vars)-1]], val = [1, -400]),
                        cplex.SparsePair(ind = [second_backward_vars[i],second_backward_vars[len(second_backward_vars)-1]], val = [1, -400])],

            senses = ["L","L"],
            rhs=[0.01,0.01],
            range_values = [0,0],
            names=["coupled_"+second_forward_vars[i], "coupled_"+second_backward_vars[i]])
        
        i += 1

    return cpx,loc_bio_1,loc_bio_2,first_forward_vars,first_backward_vars,second_forward_vars,second_backward_vars


def pFBA_joint_model(cpx,loc_bio_1,loc_bio_2,cplex,first_forward_vars,first_backward_vars,second_forward_vars,second_backward_vars):
    #Performs parsimonious Flux Balance Analysis for the joint model
    cpx.solve() #first step -> FBA 

    ## pFBA 
    bio1_val = cpx.solution.get_values(loc_bio_1)
    bio2_val = cpx.solution.get_values(loc_bio_2)

    cpx.variables.set_lower_bounds(loc_bio_1,bio1_val)
    cpx.variables.set_lower_bounds(loc_bio_2,bio2_val)

    summation = bio1_val + bio2_val
    cpx.linear_constraints.add(lin_expr = [cplex.SparsePair(ind = [first_forward_vars[len(first_forward_vars)-1],
                                                                   second_forward_vars[len(second_forward_vars)-1]],
                                                            val = [1, 1])],
                               senses = ["E"],
                               rhs = [summation],
                               range_values = [0],
                               names=["Summation_of_biomass_formation_values"])
     
    
    var_num=cpx.variables.get_num()

    i=0
    while i<var_num:
        cpx.objective.set_linear(i, -1.0)
        i+=1

    cpx.solve() #second step -> minimization of internal fluxes
    ## pFBA

    return cpx
    

    solution_status.append(cpx.solution.get_status())
    
    
    indv_sol_second = second_model_indv.optimize()

    
    gr_indv_second.append(indv_sol_second.objective_value)

    
    gr_in_pair_candida.append(cpx.solution.get_values(loc_bio_1))
    gr_in_pair_second.append(cpx.solution.get_values(loc_bio_2))

    list_of_cpx.append(cpx)

    ii += 1


def interaction_type(ccf,bcf,math):
    #Specifies the interaction type between Candida albicans and bacterial species     
    pb = math.log2(1.1)
    ib = math.log2(0.9)
    
    if ccf>pb and bcf>pb:
        int_type = 'Mutualism [+/+]'
        role = 'promoter'
        
    elif ccf<ib and bcf<ib:
        int_type = 'Competition [-/-]'
        role = 'inhibitor'
        
    elif (ccf>ib and ccf<pb) and (bcf>ib and bcf<pb):
        int_type = 'Neutralism [o/o]'
        role = 'neutral'
        
    elif (ccf>ib and ccf<pb) and bcf<ib:
        int_type = 'Amensalism [o/-]'
        role = 'neutral'
        
    elif (ccf>ib and ccf<pb) and bcf>pb:
        int_type = 'Comensalism [o/+]'
        role = 'neutral'
        
    elif ccf<ib and (bcf>ib and bcf<pb):
        int_type = 'Amensalism [-/o]'
        role = 'inhibitor'
        
    elif ccf>pb and (bcf>ib and bcf<pb):
        int_type = 'Comensalism [+/o]'
        role = 'promoter'
        
    elif ccf>pb and bcf<ib:
        int_type = 'Parasitism [+/-]'
        role = 'promoter'
        
    elif ccf<ib and bcf>pb:
        int_type = 'Parasitism [-/+]'
        role = 'inhibitor'

    return int_type,role

