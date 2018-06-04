#
# WALTerReader.py $
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

import random as rd
import os
import numpy as np
import pandas as pd


class WALTerReader():
    
    def __init__(self):
        self.fileGeno = "dict3.txt" # file with the population genotype
        self.filePheno = "data_pheno3.csv" # file with the population phenotype
        # self.dominance = "complete"
    
    def readFilePheno(self):
        file = open(self.filePheno,"rb")
        df = pd.read_csv(file,delimiter="\t")
        file.close()
        data = np.array(df)
        header = list(df)[1:np.shape(data)[1]]
        Nb = np.shape(data)[0]
        return([data,header,Nb])
        
        
    def setConversion(self,data): 
    # data is the QTL table
        N_qtl = np.shape(data)[0]
        tab = {}
         # Same dominance diet for all the QTL 
        for qtl in range(N_qtl):
            subdict = {}
            subdict["[0,0]"] = 0
        #     if self.dominance == "complete" :
        #         subdict["[0,1]"] = data[qtl,3]
        #         subdict["[1,0]"] = data[qtl,3]

       #      if self.dominance == "partial" :
        #         subdict["0,1"] == data[qtl,3]/2
        #         subdict["1,0"] == data[qtl,3]/2
        #     
        #     if self.dominance == "recessive":
        #     subdict["0,1"] == 0
        #     subdict["1,0"] == 0
        #     
        #     subdict["[1,1]"] = data[qtl,3]
            
            # Dominance diet specified for each qtl in dataQTL 
            if data[qtl,5] == "complete":
                subdict["[0,1]"] = data[qtl,3]
                subdict["[1,0]"] = data[qtl,3]
            if data[qtl,5] == "partial":
                subdict["[0,1]"] = 0.5*data[qtl,3]
                subdict["[1,0]"] = 0.5*data[qtl,3]
            if data[qtl,5] == "recessive":
                subdict["[0,1]"] = 0
                subdict["[1,0]"] = 0
            
            subdict["[1,1]"] = data[qtl,3]
            tab[qtl] = subdict
        return(tab)
        
    
    def getGenoPop(self,data):
    # data is the QTL table
        pheno = self.readFilePheno()
        N_qtl = np.shape(data)[0]
        N_ind = pheno[2]
        try :
            file = open(self.fileGeno,"rb")
            listGeno = eval(file.read())
            file.close()
            if len(list(listGeno.values())[0]) != N_qtl:
                print("Warning : incorrect file, genotype randomly generated (wrong number of qtls)")
                listGeno = {}
                for ind in range(N_ind):
                    subdict = {}
                    for qtl in range(N_qtl):
                        subdict[qtl] = [rd.choice([0,1]),rd.choice([0,1])]
                    listGeno[ind] = subdict
        except FileNotFoundError :
            listGeno = {}
            for ind in range(N_ind):
                subdict = {}
                for qtl in range(N_qtl):
                    subdict[qtl] = [rd.choice([0,1]),rd.choice([0,1])]
                listGeno[ind] = subdict
            print("Warning :incorrect file, genotype randomly generated (wrong file or path name)")
        return(listGeno)
    
    def getAllelTab(self,listGeno):
    # listGeno is the dictionary list of the population genotype
        N_ind = len(listGeno)
        N_allel = len(listGeno[0])
        tab_out = np.zeros((N_ind,N_allel),dtype=object)
        for i in range(N_ind):
            tab_out[i,:] = list(listGeno[i].values())
        return(tab_out)
    
    def getQTLtab(self,popGeno,data):
        N_ind = len(popGeno)
        conv = self.setConversion(data)
        N_qtl = len(conv)
        tab_out = np.zeros((N_ind,N_qtl),dtype=object)
        for i in range(N_ind):
            qtl_val = []
            for j in range(N_qtl):
                allel = popGeno[i][j] # value of allele i for qtl j
                # conv[j] is the list of alleles and values for qtl j
                list_allel = list(conv[j].keys()) # list of alleles for qtl j
                list_valqtl = list(conv[j].values()) # list of qtl j values for each allele
                idx = list_allel.index(allel) # allele indice in the list
                val = list_valqtl[idx] # qtl value matching this allele
                qtl_val.append(val)
            tab_out[i,:] = qtl_val
        return(tab_out)
    
    def dataPheno(self,data):
        dict_genoPop = self.getGenoPop(data)
        genoPop = self.getAllelTab(dict_genoPop)
        conv = self.setConversion(data)
        pheno = self.readFilePheno()
        return([genoPop,conv,pheno])

