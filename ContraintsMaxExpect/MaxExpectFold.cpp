//
//  MaxExpectFold.cpp
//  ContraintsMaxExpect
//
//  Created by ding xin qiang on 11/26/12.
//  Copyright (c) 2012 ding xin qiang. All rights reserved.
//

#include "MaxExpectFold.h"

// Define the constructor
MaxExpectFold::MaxExpectFold(probabilityclass &prob){
    // allocate the V, W, P array and dot
    // allocate V
    V = new double *[prob.numofbases];
    for (i = 0; i < prob.numofbases; i++){
        V[i] = new double [prob.numofbases];
    }
    
    // allocate W
    W = new double *[prob.numofbases];
    for (i = 0; i < prob.numofbases; i++){
        W[i] = new double [prob.numofbases];
    }
    
    // allocate P
    P = new int *[prob.numofbases];
    for (i = 0; i < prob.numofbases; i++){
        P[i] = new int [prob.numofbases];
    }
    
    // allocate dot
    dot = new char [prob.numofbases];
    
    // fill the V, W and P array
    for (k = 1; k <=prob.numofbases; k++){
        if (k <= 4) {
            for (i = 0; i <= prob.numofbases - k; i++) {
                V[i][i+k-1] = -1;
                P[i][i+k-1] = -1;
                W[i][i+k-1] = 0;
                for (m = i; m <= i+k-1; m++) {
                    W[i][i+k-1] += prob.sgprob[m];
                }
            }
        }else{
            for (i = 0; i <= prob.numofbases - k; i++) {
                if (! ifpair(prob.nucleotide[i], prob.nucleotide[i+k-1])){
                    V[i][i+k-1] = -1;
                }else{
                    V[i][i+k-1] = 2 * r * prob.bpprob[i][i+k-1] + W[i+1][i+k-2];
                }
                W[i][i+k-1] = V[i][i+k-1];
                P[i][i+k-1] = 0;
                for (l = i; l < i+k-1; l++) {
                    if (W[i][i+k-1] < (W[i][l] + W[l+1][i+k-1])) {
                        W[i][i+k-1] = W[i][l] + W[l+1][i+k-1];
                        P[i][i+k-1] = l;
                    }
                }
            }
        }
    }
/*    // Print the V, W and P array
	for (j = prob.numofbases-1; j >= 0; j--){
		for (i = 0; i <=j ; i++){
			cout << V[i][j] << ' ';
		}
		cout << '\n';
	}
    
	for (j = prob.numofbases-1; j >= 0; j--){
		for (i = 0; i <=j ; i++){
			cout << W[i][j] << ' ';
		}
		cout << '\n';
	}
    
	for (j = prob.numofbases-1; j >= 0; j--){
		for (i = 0; i <=j ; i++){
			cout << P[i][j] << ' ';
		}
		cout << '\n';
	}
 */
   BracketNotation(P, 0, prob.numofbases-1); 
}


// Define the function ifpair: function used to see if two nucleotides can form cannonical pairing
bool MaxExpectFold::ifpair(char baseOne, char baseTwo){
    baseOne = toupper(baseOne);
	baseTwo = toupper(baseTwo);
	switch(baseOne){
		case 'A':{
			switch (baseTwo){
				case 'A': return false;
				case 'U': return true;
				case 'G': return false;
				case 'C': return false;
			}
		}
		case 'U':{
			switch (baseTwo){
				case 'A': return true;
				case 'U': return false;
				case 'G': return true;
				case 'C': return false;
			}
		}
		case 'C':{
			switch (baseTwo){
				case 'A': return false;
				case 'U': return false;
				case 'G': return true;
				case 'C': return false;
			}
		}
		case 'G':{
			switch (baseTwo){
				case 'A': return false;
				case 'U': return true;
				case 'G': return false;
				case 'C': return true;
			}
		}
    }
}

// Define the BracketNotation function: calculate the bracket notation for RNA secondary structures according to P array and put it in dot

void MaxExpectFold::BracketNotation(int **P, int i, int j){
	if (i<j){
		if(P[i][j]==0){
			strncat(dot, "(", 1);
			BracketNotation(P,i+1,j-1);
			strncat(dot, ")", 1);
		}
		else if(P[i][j]==-1){
			strncat(dot, ".", 1);
			BracketNotation(P,i+1,j-1);
			strncat(dot, ".", 1);
		}
		else{
			BracketNotation(P,i,P[i][j]);
			BracketNotation(P,P[i][j]+1,j);
		}
	}
	else if (i==j){
		strncat(dot, ".", 1);
		return;
	}
}
