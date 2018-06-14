N_ind = 10000
N_qtl = 2

dict={}
rloc1v0=[]
rloc1v1=[]
rloc2v0=[]
rloc2v1=[]

for ind in range(N_ind):
    subdict={}
    a=rd.choice([0,1])
    b=rd.choice([0,1])
    c=rd.choice([0,1])
    d=rd.choice([0,1])
    subdict[1] = [a,b] # qtl 1
    subdict[2] = [c,d] # qtl 2
    rloc1v0.append(a)
    rloc1v1.append(b)
    rloc2v0.append(c)
    rloc2v1.append(d)
    dict[ind] = subdict
    
from scipy.stats.stats import pearsonr
loc1=rloc1v0+rloc1v1
loc2=rloc2v0+rloc2v1   
print(pearsonr(loc1,loc2))

##

res_dict_file = open("Dict_C.txt","w")
res_dict_file.write(str(dict))
res_dict_file.close()