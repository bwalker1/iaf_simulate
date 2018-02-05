//
//  spikedata.hpp
//  TE_DT
//
//  Created by Benjamin Walker on 8/24/17.
//
//

#ifndef spikedata_hpp
#define spikedata_hpp

#include <stdio.h>
#include "DTSource.h"
#include "DTRandom.h"
#include "DTRandom.h"
#include "DTIndex.h"
#include "DTDoubleArrayRegion.h"
#include "DTSource.h"
#include "DTArguments.h"
#include "DTSaveError.h"

#include "DTDoubleArray.h"
#include "DTMatlabDataFile.h"
#include "DTNetwork.h"

// Common utilities
#include "DTDoubleArrayOperators.h"
#include "DTProgress.h"
#include "DTTimer.h"
#include "DTUtilities.h"
#include "DTDictionary.h"

#include <math.h>
#include <set>


#include <queue>

typedef struct
{
    double arrival_time;
    int target;
    int source;
    double value;
} iaf_signal;

struct comp
{
    bool operator()(const iaf_signal &l, const iaf_signal &r)
    {
        return l.arrival_time > r.arrival_time;
    }
};

DTDoubleArray spikedataComputation(double n,DTDoubleArray _mu,double _alpha,double _beta,const DTDoubleArray &conn,double tmax, const DTDoubleArray &seed_list, double seed);
DTDoubleArray spikedataComputation(double n,DTDoubleArray _mu,double _alpha,double _beta,const DTDoubleArray &conn,double tmax, const DTDoubleArray &seed_list, double seed, set<int> &observed);

#endif /* spikedata_hpp */
