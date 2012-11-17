//
//  force_function.cpp
//  ContraintsMaxExpect
//
//  Created by ding xin qiang on 11/17/12.
//  Copyright (c) 2012 ding xin qiang. All rights reserved.
//

#include "force_function.h"


//force is used to enforce the folding constraints specified in a structure
//The following definitions are bitwise applied to the fce array:
//SINGLE applies to any fce(i,j) s.t. i or j should be single stranded
//PAIR applies to any fce(i,j) where i is paired to j
//NOPAIR applies to any fce(i,j) where either i or j is paired to
//	another nucleotide or i and j are forbidden to pair
//DUBLE applies to any fce(i.j) where an nuc, k, i<k<j is double stranded
//INTER applies to any fce(i,j) such that some k, i<k<j, is the virtual linker
//	used for intermolecular folding
//INTER applies to any fce(i,i) s.t. i is the center of the virtual linker
//The above terms are defined in define.h


//mod[i] is a bool array that is set to true if i is a nuc accessible to chemical
//	modification
//lfce[i] is a bool array that is set to true if i is double stranded


//force also double the sequence s.t. ct->numseq[i+N] = ct->numseq[i] where
//	N is the number of nucleotides


void force(structure *ct,forceclass *fce, bool *lfce){
	int i,j;
	register int number;
    
	number = ct->numofbases;
    
	for (i=1;i<=ct->nnopair;i++) {
		
		if(ct->nopair[i]<=ct->numofbases)	forcesingle(ct->nopair[i],ct,fce);
	}
    
	for (i=1;i<=ct->npair;i++) {
		if(ct->pair[i][0]<=ct->numofbases&&ct->pair[i][1]<=ct->numofbases) {
			forcepair(ct->pair[i][0],ct->pair[i][1],ct,fce);
			forcedbl(ct->pair[i][0],ct,fce,lfce);
			forcedbl(ct->pair[i][1],ct,fce,lfce);
		}
	}
    
	for (i=1;i<=ct->ndbl;i++) {
		if (ct->dbl[i]<=ct->numofbases) forcedbl(ct->dbl[i],ct,fce,lfce);
	}
    
    
	//u's in gu pairs must be double stranded
	for (i=0;i<ct->ngu;i++) {
		if (ct->gu[i]<=ct->numofbases) forcedbl(ct->gu[i],ct,fce,lfce);
	}
    
    
    
	if (ct->intermolecular) {//this indicates an intermolecular folding
		//for (i=0;i<=3;i++) {
        // 	forcesingle(ct->inter[i])//don't allow the intermolecular indicators to pair
        //}
        for (i=0;i<3;i++) {
            forceinter(ct->inter[i],ct,fce);
        }
        
        
        
        
        fce->f(ct->inter[1],ct->inter[1]) = fce->f(ct->inter[1],ct->inter[1])|INTER;
        
        
	}
    
	for (i=0;i<ct->nforbid;i++) {
		if(ct->forbid[i][0]<=ct->numofbases&&ct->forbid[i][1]<=ct->numofbases) fce->f(ct->forbid[i][0],ct->forbid[i][1]) = fce->f(ct->forbid[i][0],ct->forbid[i][1])|NOPAIR;
		if(ct->forbid[i][0]<=ct->numofbases&&ct->forbid[i][1]<=ct->numofbases) fce->f(ct->forbid[i][1],ct->forbid[i][0]+ct->numofbases)=fce->f(ct->forbid[i][1],ct->forbid[i][0]+ct->numofbases)|NOPAIR;
	}
    
	//Double up the sequence
	for (i=1;i<=number;i++) {
        ct->numseq[(number)+i] = ct->numseq[i];
	}
    
	//The next section handles the case where base pairs are not
    //not allowed to form between nucs more distant
    //than ct->maxdistance
	if (ct->limitdistance) {
        
		if (!ct->templated) ct->allocatetem();
        
		for (j=minloop+2;j<=ct->numofbases;j++) {
			for (i=1;i<j;i++) {
				if (j-i>=ct->maxdistance) ct->tem[j][i]=false;
			}
		}
        
        
		
	}
    
    
}

//this function is used to set up the fce array for a base, x, that should be single-stranded
void forcesingle(int x,structure* ct,forceclass *v) {
	int i;
    
    for (i=x;i<x+(ct->numofbases);i++) {
        v->f(x,i)=v->f(x,i)|SINGLE;
    }
    for (i=1;i<=x;i++) {
        v->f(i,x)=v->f(i,x)|SINGLE;
    }
    for (i=x+1;i<=ct->numofbases;i++) {
        v->f(i,x+ct->numofbases) = v->f(i,x+ct->numofbases)|SINGLE;
    }
}





void forcepair(int x,int y,structure *ct,forceclass *v) {
	int i,j;
    v->f(x,y) = v->f(x,y)|PAIR;
    v->f(y,x+ct->numofbases)=v->f(y,x+ct->numofbases)|PAIR;
    for (i=y+1;i<=x-1+ct->numofbases;i++) {
        v->f(x,i) = v->f(x,i)|NOPAIR;
    }
    for (i=x;i<=y-1;i++) {
        v->f(x,i) = v->f(x,i)|NOPAIR;
    }
    for (i=1;i<=x-1;i++) {
        v->f(i,y) = v->f(i,y)|NOPAIR;
    }
    for (i=x+1;i<=y;i++) {
        v->f(i,y) = v->f(i,y)|NOPAIR;
    }
    for (i=1;i<=x-1;i++) {
        v->f(i,x) = v->f(i,x)|NOPAIR;
    }
    for (i=y+1;i<=ct->numofbases;i++) {
        v->f(i,y+ct->numofbases)=v->f(i,y+ct->numofbases)|NOPAIR;
    }
    for (i=y;i<=x-1+(ct->numofbases);i++) {
        v->f(y,i) = v->f(y,i)|NOPAIR;
    }
    for (i=(ct->numofbases)+x+1;i<=(ct->numofbases)+y-1;i++) {
        v->f(y,i) = v->f(y,i)|NOPAIR;
    }
    for (i=x+1;i<=y-1;i++) {
        v->f(i,x+ct->numofbases) = v->f(i,x+ct->numofbases)|NOPAIR;
    }
    for (i=y+1;i<=ct->numofbases;i++) {
        v->f(i,x+ct->numofbases) = v->f(i,x+ct->numofbases)|NOPAIR;
    }
    for (i=1;i<=x-1;i++) {
        for (j = x+1;j<=y-1;j++){
            v->f(i,j) = v->f(i,j)|NOPAIR;
        }
    }
    for (i=x+1;i<=y-1;i++) {
        for (j=y+1;j<=(ct->numofbases)+x-1;j++) {
            v->f(i,j) = v->f(i,j)|NOPAIR;
        }
    }
    for (i=y+1;i<=ct->numofbases;i++) {
        for (j=(ct->numofbases)+x+1;j<=(ct->numofbases)+y-1;j++) {
            v->f(i,j) = v->f(i,j)|NOPAIR;
        }
    }
}



void forcedbl(int dbl,structure* ct,forceclass *w,bool *v) {
	int i,j;
    
	
    v[dbl] = true;
    v[dbl+ct->numofbases] = true;
    
    
    for(i=dbl+1;i<=ct->numofbases;i++) {
        for (j=1;j<dbl;j++) {
            w->f(j,i) = w->f(j,i)|DUBLE;
        }
    }
    for(j=(dbl+(ct->numofbases)-1);j>ct->numofbases;j--) {
        for (i=dbl+1;i<=ct->numofbases;i++) {
            w->f(i,j) = w->f(i,j)|DUBLE;
        }
    }
    
    
}

void forceinter(int dbl,structure* ct,forceclass *w) {
	int i,j;
    
    
    for(i=dbl+1;i<=ct->numofbases;i++) {
        for (j=1;j<dbl;j++) {
            w->f(j,i) = w->f(j,i)|INTER;
        }
    }
    for(j=(dbl+(ct->numofbases)-1);j>ct->numofbases;j--) {
        for (i=dbl+1;i<=ct->numofbases;i++) {
            w->f(i,j) = w->f(i,j)|INTER;
        }
    }
    for(i=dbl+1+ct->numofbases;i<=2*ct->numofbases;i++) {
        for (j=ct->numofbases;j<dbl+ct->numofbases;j++) {
            w->f(j,i) = w->f(j,i)|INTER;
        }
    }
    
    
    
    
    
}

void forceinterefn(int dbl,structure* ct,forceclass *w) {
	int i,j;
    
    
    for(i=dbl+1;i<=ct->numofbases;i++) {
        for (j=1;j<dbl;j++) {
            w->f(j,i) = w->f(j,i)|INTER;
        }
    }
}

//returns true if pair i and j is not a GU pair and the most dajacent pairs are not GU
//	used with chemical modification data to make sure that some contributions are not over-counted
bool notgu(int i, int j, structure *ct) {
    
	if (ct->numseq[i]==3&&ct->numseq[j]==4) return false;
	else if (ct->numseq[i]==4&&ct->numseq[j]==3) return false;
	else if (ct->numseq[i+1]==3&&ct->numseq[j-1]==4) return false;
	else if (ct->numseq[i+1]==4&&ct->numseq[j-1]==3) return false;
	else if (ct->numseq[i-1]==3&&ct->numseq[j+1]==4) return false;
	else if (ct->numseq[i-1]==4&&ct->numseq[j+1]==3) return false;
	else return true;
    
}