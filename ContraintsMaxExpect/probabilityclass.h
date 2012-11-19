//
//  probabilityclass.h
//  ContraintsMaxExpect
//
//  Created by ding xin qiang on 11/19/12.
//  Copyright (c) 2012 ding xin qiang. All rights reserved.
//

#ifndef __ContraintsMaxExpect__probabilityclass__
#define __ContraintsMaxExpect__probabilityclass__

/*
 The probabilityclass encapulates the probability of basepairing and single strandness for a RNA secondary structure. The class is constructed using a .pfs file.
*/
#include <iostream>

#include "structure.h"
#include "defines.h"
#include "pfunction.h"

class probabilityclass{
    public:
    // numofbases holds the num of nucleotide in the sequence
    int numofbases;
    // bpprob[][] holds the basepairing probability of the structure.
    double **bpprob;
    // sgprob[] holds the single-strandness probability of each nucleotide.
    double *sgprob;
    // nucleotide[] holds the nucleotide of the sequence
    char *nucleotide;
    // Constructor to build the class using .pfs file name
    probabilityclass(char *pfs_filename);
   
};
#endif /* defined(__ContraintsMaxExpect__probabilityclass__) */
