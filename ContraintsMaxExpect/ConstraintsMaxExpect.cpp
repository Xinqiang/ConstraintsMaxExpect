//
//  main.cpp
//  ContraintsMaxExpect
//
//  Created by ding xin qiang on 11/17/12.
//  Copyright (c) 2012 ding xin qiang. All rights reserved.
//

#include <iostream>
#include <cstring>
#include "probabilityclass.h"
#include "MaxExpectFold.h"

using namespace std;

int main(int argc, char * argv[])
{
    // Build the probabilitycalss prob using argv[1]: .pfs file name
    probabilityclass prob(argv[1]);
    
    // build the MaxExpectFold class maxexpect using prob
    MaxExpectFold maxexpect(prob);
    
    
    cout << prob.nucleotide << "\n";
    cout << maxexpect.dot;
}