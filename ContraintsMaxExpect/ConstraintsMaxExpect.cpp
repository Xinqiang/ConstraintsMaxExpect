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
    structure *ct;
    double *w5;
    double *w3;
    pfunctionclass *v;
    pfunctionclass *w;
    pfunctionclass *wmb;
    pfunctionclass *wl;
    pfunctionclass *wmbl;
    pfunctionclass *wcoax;
    forceclass *fce;
    double *scaling;
    bool *mod;
    bool *lfce;
    pfdatatable *data;
    
    readpfsave("/Users/xinqiang/Research/RNASecondaryStructure/ContraintsMaxExpect/ContraintsMaxExpect/test.pfs", ct, w5, w3, v, w, wmb, wl, wmbl, wcoax, fce, scaling, mod, lfce, data);
    
    cout << w3[1];
    cout << "Hello, World!\n";
    return 0;
}

