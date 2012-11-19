//
//  probabilityclass.cpp
//  ContraintsMaxExpect
//
//  Created by ding xin qiang on 11/19/12.
//  Copyright (c) 2012 ding xin qiang. All rights reserved.
//

#include "probabilityclass.h"
using namespace std;

probabilityclass::probabilityclass(char *pfs_filename){
    
    // Open the .pfs file
    ifstream pfsin(pfs_filename, ios::in | ios::binary);
    
    // Read the .pfs file
    if(!pfsin){
        cout << "Cannot open file.\n";
    }else{
        
        // read the version and numofbase in order to initiallize the structure
        short vers;
        //int  numofbases;
        pfsin.read((char *) &vers, sizeof vers);
        structure *ct = new structure;
        pfsin.read((char *) &ct->numofbases, sizeof numofbases);
        
        numofbases = ct->numofbases;
        // Initiallize the structure ct using numofbase.
        ct->allocate(ct->numofbases);
        
        // Export the nucleotides to the public variable of the class probabilitycalss.
        nucleotide = &(ct->nucs[1]);
        
        // Initiallize the array to holds the information from .pfs file.
        double *w5 = new double [ct->numofbases + 1];
        double *w3 = new double [ct->numofbases + 2];
        pfunctionclass *v = new pfunctionclass(ct->numofbases);
        pfunctionclass *w = new pfunctionclass(ct->numofbases);
        pfunctionclass *wmb = new pfunctionclass(ct->numofbases);
        pfunctionclass *wl = new pfunctionclass(ct->numofbases);
        pfunctionclass *wmbl = new pfunctionclass(ct->numofbases);
        pfunctionclass *wcoax = new pfunctionclass(ct->numofbases);
        forceclass *fce = new forceclass(ct->numofbases);
        double *scaling = new double;
        bool *mod = new bool [2*ct->numofbases + 1];
        bool *lfce = new bool [2*ct->numofbases + 1];
        pfdatatable *data = new pfdatatable();
        
        // Now read the information from .pfs file
        readpfsave(pfs_filename, ct, w5, w3, v, w, wmb, wl, wmbl, wcoax, fce, scaling, mod, lfce, data);
        
    // Using the information from .pfs file to calculate the pairing probability using calculate probability
        // Initiallize the array[][] to holds the base pairing probability
        // Note! The index of the array is not the same with the the position of nucleotide in the sequence
        bpprob = new double *[ct->numofbases];
        int k;
        for (k=0;k<ct->numofbases;k++){
            bpprob[k] = new double [ct->numofbases];
        }
        
        //Initaillize the array [] to holds the single strandness probability
        sgprob = new double [ct->numofbases];
        
        // Calculate the base pairing probability and single strandness probability using calculateprobability() and store them in bppro[][] and sgprob[]
        // i and j are the index of position of nucleotide in the sequence
        int i,j;
        for (i=0; i< ct->numofbases; i++){
            for (j=0; j< ct->numofbases; j++){
                bpprob[i][j] = calculateprobability(i+1,j+1, v, w5, ct, data, lfce, mod, *scaling, fce);
            }
        }
        
        // Calculate the single strandness of probability
        // sumofbsprob stores the sum of base pairing probability for a nucleotide
        double sumofbsprob;
        for (i = 0; i < ct->numofbases; i++){
            sumofbsprob = 0;
            for (j = i + 1; j < ct->numofbases; j++){
                sumofbsprob += bpprob[i][j];
            }
            for (k = 0; k < i; k++){
                sumofbsprob += bpprob[k][i];
            }
            sgprob[i] = 1 - sumofbsprob;
        }
    }
}