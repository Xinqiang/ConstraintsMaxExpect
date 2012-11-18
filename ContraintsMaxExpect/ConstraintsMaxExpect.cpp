//
//  main.cpp
//  ContraintsMaxExpect
//
//  Created by ding xin qiang on 11/17/12.
//  Copyright (c) 2012 ding xin qiang. All rights reserved.
//

#include <iostream>
#include "structure.h"
#include "defines.h"
#include "pfunction.h"

using namespace std;
int main(int argc, const char * argv[])
{
    ifstream pfsin("/Users/xinqiang/Research/RNASecondaryStructure/ContraintsMaxExpect/ContraintsMaxExpect/test.pfs", ios::in | ios::binary);
    
    if(!pfsin){
        cout << "Cannot open file.\n";
        return 1;
    }
    short vers;
    int  numofbases;
    
    pfsin.read((char *) &vers, sizeof vers);
    
    cout << vers << "\n"; 
    
    structure *ct = new structure;
    
    pfsin.read((char *) &ct->numofbases, sizeof numofbases);
    
    cout << ct->numofbases << "\n";
    
    ct->allocate(ct->numofbases);
    
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
    
    readpfsave("/Users/xinqiang/Research/RNASecondaryStructure/ContraintsMaxExpect/ContraintsMaxExpect/test.pfs", ct, w5, w3, v, w, wmb, wl, wmbl, wcoax, fce, scaling, mod, lfce, data);
    
    cout << w3[1] << "\n" << ct->nucs << "\n";
    cout << *scaling << "\n";
    int i,j;
    
    for (i=1; i<=ct->numofbases; i++){
        for (j=i+1; j<=ct->numofbases; j++){
            cout << calculateprobability(i,j, v, w5, ct, data, lfce, mod, *scaling, fce) << " ";
        }
        cout << "\n";
    }
    return 0;
}

