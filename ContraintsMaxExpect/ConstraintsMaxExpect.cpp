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

using namespace std;

int main(int argc, char * argv[])
{
    // Build the probabilitycalss prob using argv[1]: .pfs file name
    probabilityclass prob(argv[1]);
    
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
    
    // V[i][j] holds the Expected Accuracy of the segment between i and j on the condition i and j form a base pair.
    //  V[i][j] = -1, if j-i+1 <=4
    //  V[i][j] = -1, if i and j cannot form cannonical pair.
    double **V;
    V = new double *[prob.numofbases];
    for (i = 0; i < prob.numofbases; i++){
        V[i] = new double [prob.numofbases];
    }
    
    // W[i][j] holds the Expected Accuracy of the segment between i and j
    double **W;
    W = new double *[prob.numofbases];
    for (i = 0; i < prob.numofbases; i++){
        W[i] = new double [prob.numofbases];
    }
    
    // P(i,j) is the k so that W(i,j) = W(i,k) + W(k+1,j) and i <= k < j,
    // but if W(i,j) = V(i,j), then P(i,j) = 0; if j-i+1<=4, then P(i,j)=-1.
    int **P = new int *[prob.numofbases];
    for (i = 0; i < prob.numofbases; i++){
        P[i] = new int [prob.numofbases];
    }
    
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
    
    
// Print the V, W and P array
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
    
    cout << prob.nucleotide << "\n";
    BracketNotation(P, 0, prob.numofbases-1);
    return 0;
}


// Print the bracket notation for RNA secondary structures
void BracketNotation(int **P, int i, int j){
	if (i<j){
		if(P[i][j]==0){
			cout << '(';
			BracketNotation(P,i+1,j-1);
			cout << ')';
		}
		else if(P[i][j]==-1){
			cout << '.';
			BracketNotation(P,i+1,j-1);
			cout << '.';
		}
		else{
			BracketNotation(P,i,P[i][j]);
			BracketNotation(P,P[i][j]+1,j);
		}
	}
	else if (i==j){
		cout << '.';
		return;
	}
}


// function used to see if two nucleotides can form cannonical pairing
bool ifpair(char baseOne, char baseTwo){
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