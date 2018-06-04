#
# Offspring.py $
#
# This file is part of GenoWALT, a sofware of simulation of crop population.
#
# Copyright (C) 2017 Meije Gawinowski (meije.gawinowski@wanadoo.fr)
#
# 
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#




import simuPOP as sim
import numpy as np
import csv
import copy as cp





class Offspring():


    def __init__(self,trt_interest):
        self.trt_interest = trt_interest
        self.trt_mean = [51.9]
        self.choice_pheno = "rescaling"


    def offspringGeno(self,offpop,conv,tab_qtl):
        """ This function returns the allel tab for the offspring population """
        # indices des QTL
        qtl_names = tab_qtl[:,4]
        qtl_indices = list(offpop.lociByNames(qtl_names))
        
        N_ind = offpop.popSize()
        N_qtl = len(qtl_indices)
        
        tab_out = np.zeros((N_ind,N_qtl),dtype=object)
        for i in range(N_qtl):
            list_allel = []
            for j in range(N_ind):
                ind = offpop.individual(j)
                val = [int(ind.allele(qtl_indices[i],0)),int(ind.allele(qtl_indices[i],1))] # for each individual j we take the QTL i and we take the allele value on the homoeologous copy 0 of the chromosome 
                list_allel.append(val)
            tab_out[:,i] = list_allel
        print("genotab",tab_out)
        return(tab_out)
        
    def dictGenoOffspring(self,genotab):
        dict = {}
        for i in range(np.shape(genotab)[0]):
            subdict = {}
            for j in range(np.shape(genotab)[1]):
                subdict[j] = str(genotab[i,j])
            dict[i] = subdict
        print("offspring geno dict",dict)
        return(dict)

    def offspringQTLTab(self,dictGeno,conv):
        N_ind = len(dictGeno)
        N_qtl = len(conv)
        tab_out = np.zeros((N_ind,N_qtl),dtype=object)
        for i in range(N_ind):
            qtl_val = []
            for j in range(N_qtl):
                allel = dictGeno[i][j] # allel value for the individual i at the qtl j       
                al2="["+str(allel[1])+","+str(allel[4])+"]"
                # conv[j] est liste allèles/valeurs pour le qtl j
                list_allel = list(conv[j].keys()) # alleles list for qtl j
                list_valqtl = list(conv[j].values()) # qtl j values list for each allele
                try :
                    idx = list_allel.index(allel) # index of the allele in the list
                except ValueError :
                    idx = list_allel.index(al2)
                val = list_valqtl[idx] # qtl value matching this allele
                qtl_val.append(val)
            tab_out[i,:] = qtl_val
        return(tab_out)
    
    
    def rescaling(self,init,real_min,real_max):
        init_min = min(init)
        if init_min < 0:
            init = list(np.add(-init_min,init))
        y = [(real_max-real_min)/(max(init)-min(init))*init[i]+real_min for i in range(len(init))]
        return(y)
    
    def trt_indices2(self,tab_qtl,trt):
            sublist = []
            for i in range(np.shape(tab_qtl)[0]):
                list_trt = tab_qtl[i,0]
                try :
                    idx = list_trt.index(trt)
                    sublist.append(idx)
                except ValueError :
                    sublist.append(-1)
            return(sublist)
    
    
    def offspringPheno(self,offQTLtab,choice,tab_qtl):
        N_ind = np.shape(offQTLtab)[0]
        N_qtl = np.shape(offQTLtab)[1]
        tab_pheno = np.zeros((N_ind,1),dtype=object)
        tab_pheno[:,0] = list(range(N_ind))
        header = ['individual'] + self.trt_interest
        print("TABLE",offQTLtab)
        for trt in self.trt_interest :
            trt_idx = self.trt_interest.index(trt)
            list_trt_indices = self.trt_indices2(tab_qtl,trt)
            print("list",list_trt_indices)
            newQTLtab = cp.deepcopy(offQTLtab)
            comp = 0
            qtl = 0
            while comp <= N_qtl-1:
                # print("qtl",qtl)
                if list_trt_indices[comp] == -1:
                    newQTLtab = np.delete(newQTLtab,qtl,1)
                    comp = comp + 1
                else:
                    for ind in range(N_ind):
                        # print("ind",ind)
                        # print("shape",np.shape(newQTLtab))
                        list_valqtl = newQTLtab[ind,qtl]
                        valqtl_idx = list_trt_indices[qtl]
                        # print("idx",valqtl_idx)
                        # print(newQTLtab[ind,qtl])
                        newQTLtab[ind,qtl] = newQTLtab[ind,qtl][valqtl_idx]
                    qtl = qtl + 1
                    comp = comp + 1
                    
            if choice == "file" :
                val_trt = [self.trt_mean[trt_idx]+sum(newQTLtab[i,:]) for i in range(N_ind)]
            if choice == "sampling" :
                val_trt = [sum(newQTLtab[i,:]) for i in range(N_ind)]
                val_trt = self.rescaling(val_trt,57.9,121.6)
            buf_trt = np.zeros((N_ind,1),dtype=object)
            buf_trt[:,0] = val_trt
            tab_pheno = np.concatenate((tab_pheno,buf_trt),1)
            
        header_buf = np.zeros((1,len(header)),dtype=object)
        header_buf[0,:] = header
        np.concatenate((header_buf,tab_pheno),0)
        
        return(tab_pheno)
        
        
        
    def Result(self,offpop,conv,tab_qtl,choice):
        genotab = self.offspringGeno(offpop,conv,tab_qtl)
        geno_dict = self.dictGenoOffspring(genotab)
        offQTLtab = self.offspringQTLTab(geno_dict,conv)
        tabpheno = self.offspringPheno(offQTLtab,choice,tab_qtl)
        
        # Writing of a genotype dictionary in a txt file
        res_dict_file = open("dictRes.txt","w")
        res_dict_file.write(str(geno_dict))
        res_dict_file.close()
        
        # Writing of a PH phenotype table in a csv file
        res_pheno_file = open("phenoRes.csv","w")
        pheno_writer = csv.writer(res_pheno_file,delimiter="\t")
        pheno_writer.writerows(tabpheno)
        res_pheno_file.close()
       
        return([geno_dict,tabpheno])
        

    

