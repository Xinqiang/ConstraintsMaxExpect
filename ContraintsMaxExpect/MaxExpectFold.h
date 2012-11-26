//
//  MaxExpectFold.h
//  ContraintsMaxExpect
//
//  Created by ding xin qiang on 11/26/12.
//  Copyright (c) 2012 ding xin qiang. All rights reserved.
//

#ifndef __ContraintsMaxExpect__MaxExpectFold__
#define __ContraintsMaxExpect__MaxExpectFold__

// The class MaxExpectFold is used for fold RNA sequence into secondary structure using MaxExpect rules. It require the class probabilityclass as parameter to construct

#include <iostream>
#include "probabilityclass.h"
#include <string.h>

class MaxExpectFold {
    // the firt base index of a RNA segment.
    int i;
    // the end base index of a RNA segment.
	int j;
    // the base index between i and j, so i<=l<j.
	int l;
    // the length of RNA between base i and base j, so k = j-i+1.
	int k;
    // index between i and i+k-1
    int m;
    // the coefficent in the definatioin of Expected accuracy: S = 2*r*bsprob() + sgprob()
    int r = 1;
    // funcition to see if two nucleotides can form cannocial pair
    bool ifpair(char baseOne, char baseTwo);
    // Function used to generate the BracketNotation of RNA secondary structure recursively
    void BracketNotation(int **P, int i, int j);
    
public:
    // V[i][j] holds the Expected Accuracy of the segment between i and j on the condition i and j form a base pair.
    //  V[i][j] = -1, if j-i+1 <=4
    //  V[i][j] = -1, if i and j cannot form cannonical pair.
    double **V;
    
    // W[i][j] holds the Expected Accuracy of the segment between i and j
    double **W;
    
    // P(i,j) is the k so that W(i,j) = W(i,k) + W(k+1,j) and i <= k < j,
    // but if W(i,j) = V(i,j), then P(i,j) = 0; if j-i+1<=4, then P(i,j)=-1.
    int **P;
    
    // dot is a string to hold the bracket notation of the folded RNA secondary structure
    char *dot;
    
    // Constructor used to fill the V, W, P array and fold the RNA into dot. It require probabilityclass as the parameter
    MaxExpectFold(probabilityclass &prob);
    
};


#endif /* defined(__ContraintsMaxExpect__MaxExpectFold__) */



