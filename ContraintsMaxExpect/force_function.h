//
//  force_function.h
//  ContraintsMaxExpect
//
//  Created by ding xin qiang on 11/17/12.
//  Copyright (c) 2012 ding xin qiang. All rights reserved.
//

#ifndef __ContraintsMaxExpect__force_function__
#define __ContraintsMaxExpect__force_function__

#include <iostream>
#include "structure.h"
#include "forceclass.h"

//the force... functions are used to initialize the arrays used to apply
//constraints during the folding process
void forcepair(int x,int y,structure *ct,forceclass *v);
void forcesingle(int x,structure* ct,forceclass *v);
void forcedbl(int dbl,structure* ct,forceclass *w,bool *v);
void forceinter(int dbl,structure* ct,forceclass *w);
void forceinterefn(int dbl,structure* ct,forceclass *w);

//force is used to prepare arrays for the function dynamic, used during the
//	fill routines - it coordinates the force...() functions above
void force(structure *ct,forceclass *fce, bool *lfce);


bool notgu(int i, int j, structure *ct);

#endif /* defined(__ContraintsMaxExpect__force_function__) */
