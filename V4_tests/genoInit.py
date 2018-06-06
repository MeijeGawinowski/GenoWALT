N_ind = 100
N_qtl = 2
dict={}
for ind in range(N_ind):
    subdict={}
    for qtl in range(N_qtl):
        subdict[qtl] = [rd.choice([0,1]),rd.choice([0,1])]
    dict[ind] = subdict
    
res_dict_file = open("Dict0.txt","w")
res_dict_file.write(str(dict))
res_dict_file.close()