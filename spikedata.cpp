//
//  spikedata.cpp
//  TE_DT
//
//  Created by Benjamin Walker on 8/24/17.
//
//

#include "spikedata.h"
#include <assert.h>

DTDoubleArray spikedataComputation(double n,DTDoubleArray _mu,double _alpha,double _beta,const DTDoubleArray &conn,double tmax, const DTDoubleArray &seed_list, double seed)
{
    set<int> observed;
    return spikedataComputation(n, _mu, _alpha, _beta, conn, tmax, seed_list, seed, observed);
}

DTDoubleArray spikedataComputation(double n,DTDoubleArray _mu,double _alpha,double _beta,const DTDoubleArray &conn,double tmax, const DTDoubleArray &seed_list, double seed, set<int> &observed)
{
    printf("Initializing\n");

    vector<vector<int> > conn_c(n);
    for (int i=0;i<n;++i)
    {
        conn_c[i]=vector<int>();
        for (int j=0;j<n;++j)
        {
            if (conn(0,i,j)!=0)
            {
                conn_c[i].push_back(j);
            }
        }
    }
    
    
    // move forward in time exactly tmax units of time
    double t = 0;   // simulation time in current function call
    const double dt_base = 0.0001;
    
    const double threshold = 1;
    const double reset = 0;
    const bool recording = true;
    double timeLastRecord = -1000.0;
    double recordInterval = 0.001;
    
    int records = 0;
    int max_records = ceil(tmax/recordInterval);
    
    int iters = 0;
    int stride = 10;
    
    DTRandom rng(int(1000*seed_list(2)));
    
    priority_queue<iaf_signal, deque<iaf_signal>, comp> pq;
    
    DTMutableDoubleArray states(n);
    DTMutableDoubleArray mu;
    DTMutableDoubleArray alpha(n);
    DTMutableDoubleArray beta(n);

    if (_mu.Length()==n)
    {
        mu = _mu.Copy();
    }
    else if (_mu.Length()>=1)
    {
        mu = DTMutableDoubleArray(n);
        mu = _mu(0);
    }
    else
    {
        mu = DTMutableDoubleArray(n);
        mu = 1;
    }
    for (int i=0;i<n;++i)
    {
        states(i) = 0;
        alpha(i) = _alpha;
        beta(i) = _beta;
    }
    
    DTMutableDoubleArray fireFlag(n);
    DTMutableDoubleArray recordings(n,max_records);
    for (int i=0;i<n;++i)
    {
        fireFlag(i)=0;
        for (int j=0;j<max_records;++j)
        {
            recordings(i,j)=0;
        }
    }
    
    
    printf("Beginning time simulation\n");
    for (;;)
    {
        double dt;
        if (pq.empty())
        {
            dt = dt_base;
        }
        else
        {
            // simulate up to base time step, but stop early if an impulse will arrive sooner or if we're done sooner
            dt = min(min(dt_base, pq.top().arrival_time - t),tmax-t);
            dt = dt_base;
        }
        double sqrtdt = sqrt(dt);
        t += dt;
        
        // go through each neuron, applying the Integrate-and-Fire SDE through the elapsed time.
        states += _alpha*(mu-states)*dt;
        for (int i = 0; i < n; ++i)
        {
            states(i) += beta(i)*sqrtdt*rng.Normal(0.0,1.0);
        }
        
        // apply any iaf_signals that have reached their destination
        if (!pq.empty())
        {
            while (t >= pq.top().arrival_time)
            {
                int i = pq.top().target;
                states(i) += pq.top().value;
                
                
                pq.pop();
                if (pq.empty())
                {
                    break;
                }
            }
        }
        
        
        // check for firing
        for (int i = 0; i < n; ++i)
        {
            if (states(i) >= threshold)
            {
                // reset
                states(i) = reset;
                fireFlag(i) = 1;
                
                // create new outgoing iaf_signals
                for (int k = 0; k < conn_c[i].size(); ++k)
                {
                    iaf_signal out;
                    out.target = conn_c[i][k];
                    out.source = i;
                    out.value = conn(0,i,out.target);
                    out.arrival_time = t+conn(1,i,out.target);
                    pq.push(out);
                }
            }
        }
        
        
        if (recording && iters%stride==0)
        {
            // record whether each neuron has spiked recently
            recordings(DTIndex::All,records)=fireFlag;
            
            for (int i=0;i<n;++i)
            {
                fireFlag(i)=0;
            }
            timeLastRecord = t;
            
            ++records;
            if (records >= max_records)
            {
                break;
            }
        }
        ++iters;
    }
    
    printf("Time simulation complete\n");
    DTRandom rng2(seed);
    set<int> regAObserved;
    set<int> regBObserved;
    
    assert(int(n)%2 == 0);
    int cutoff = n/2;
    
    if (n>50)
    {
        while (regAObserved.size() < 25)
        {
            int val = floor(rng2.UniformHalf()*double(n)/2);
            assert(val < cutoff);
            regAObserved.insert(floor(rng2.UniformHalf()*double(n)/2));
        }
        while (regBObserved.size() < 25)
        {
            int val = int(n/2)+floor(rng2.UniformHalf()*double(n)/2);
            assert (val >= cutoff);
            regBObserved.insert(val);
        }
    }
    observed.insert(regAObserved.begin(),regAObserved.end());
    observed.insert(regBObserved.begin(),regBObserved.end());
    if (n>50)
    {
        assert(observed.size() == 50);
    }
    int nObserved;
    if (n <= 50)
    {
        nObserved = n;
        
        // for convenience
        for (int i=0;i<n;++i)
        {
            observed.insert(i);
        }
    }
    else
    {
        nObserved = observed.size();
    }
    

    
    // compress recordings
    vector<pair<int, vector<int> > > temp;
    int fires = 0;
    int observed_count = 0;
    for (int i=0;i<recordings.m();++i)
    {
        if (n>50 && !observed.count(i)) continue;
        vector<int> v1;
        for (int j=0;j<recordings.n();++j)
        {
            if (recordings(i,j)!=0)
            {
                v1.push_back(j);
                ++fires;
            }
        }
        temp.push_back(pair<int,vector<int> >(observed_count++,v1));
    }
    
    DTMutableDoubleArray toReturn(2*nObserved + fires);
    int c=-1;
    for (int i=0;i<temp.size();++i)
    {
        toReturn(++c)=temp[i].first;
        toReturn(++c)=temp[i].second.size();
        for (int j=0;j<temp[i].second.size();++j)
        {
            toReturn(++c)=temp[i].second[j];
        }
    }
    printf("Spike data complete\n");
    return toReturn;
}
