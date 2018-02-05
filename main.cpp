// Include DTSource.h if you want to include all the headers.

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

using namespace std;



#include "DTRandom.h"
#include "DTIndex.h"
#include "DTDoubleArrayRegion.h"
#include "DTSource.h"
#include "spikedata.h"
DTDoubleArray Computation(double n,double mu,double alpha,double beta,const DTDoubleArray &conn,
                          double tmax,const DTDoubleArray &seed_list,double seed,
                          const DTDoubleArray &mu_list);

int main(int argc,const char *argv[])
{
    DTSetArguments(argc,argv);
    
    DTMatlabDataFile inputFile("Input.mat",DTFile::ReadOnly);
    // Read in the input variables.
    double n = inputFile.ReadNumber("n");
    double mu = inputFile.ReadNumber("mu");
    double alpha = inputFile.ReadNumber("alpha");
    double beta = inputFile.ReadNumber("beta");
    DTDoubleArray conn = inputFile.ReadDoubleArray("conn");
    double tmax = inputFile.ReadNumber("tmax");
    DTDoubleArray seed_list = inputFile.ReadDoubleArray("seed list");
    double seed = inputFile.ReadNumber("seed");
    DTDoubleArray mu_list = inputFile.ReadDoubleArray("mu_list");
    
    // The computation.
    DTDoubleArray computed;
    clock_t t_before = clock();
    computed = Computation(n,mu,alpha,beta,conn,tmax,seed_list,seed,mu_list);
    clock_t t_after = clock();
    double exec_time = double(t_after-t_before)/double(CLOCKS_PER_SEC);
    
    // Write the output.
    DTMatlabDataFile outputFile("Output.mat",DTFile::NewReadWrite);
    
    // Output from computation
    outputFile.Save(computed,"Var");
    outputFile.Save("Array","Seq_Var");
    
    // The execution time.
    outputFile.Save(exec_time,"ExecutionTime");
    outputFile.Save("Real Number","Seq_ExecutionTime");
    
    // The errors.
    DTSaveError(outputFile,"ExecutionErrors");
    outputFile.Save("StringList","Seq_ExecutionErrors");
    
    return 0;
}

DTDoubleArray Computation(double n,double mu,double alpha,double beta,const DTDoubleArray &conn,
                          double tmax,const DTDoubleArray &seed_list,double seed,
                          const DTDoubleArray &mu_list)
{
    return spikedataComputation(n,mu_list,  alpha,  beta,   conn,  tmax,  seed_list,  seed);
}
