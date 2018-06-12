N_ind = 100000
N_qtl = 1
dict={}
for ind in range(N_ind):
    subdict={}
    for qtl in range(N_qtl):
        subdict[qtl] = [rd.choice([0,1]),rd.choice([0,1])]
    dict[ind] = subdict
    
res_dict_file = open("Dict4_1.txt","w")
res_dict_file.write(str(dict))
res_dict_file.close()