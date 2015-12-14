#ifndef FFTTOOLS_AVERAGER_H
#define FFTTOOLS_AVERAGER_H
 
#include <iostream>
#include <vector>
/* Class used to facilitate averaging of RF Waveforms 
 * Correlation is used to find the best shift 
 *
 **/ 

class TGraph; 
class TTree; 
class TEventList;

namespace FFTtools
{
    class Averager
    {


      public: 

        Averager(); 
        ~Averager() { reset(); }


        /*** Set the time window we care about (i.e. crop waveforms before averaging) **/ 
        void setTimeWindow(double tmin, double tmax){ t0 = tmin; t1 = tmax; }
        void unsetTimeWindow(){ t0 = 0; t1= 0;} 


        /*** Set the maximum amount of shift allowed **/ 
        void setMaxShift(double ms) { max_shift = ms; } 

        /*** Add g to average; assumed g is interpolated and same size as first entry  */ 
        void add(const TGraph *g); 

        /*** Interpolate g and then add to average. If interpolate_dt is 0, same as above. 
         **/
        void add(const TGraph *g, double interpolate_dt); 

        /*** Add all entries from a tree branch in event list. if list is NULL, all entries will be added.  Graphs are interpolated according to interpolation_dt  */ 
        void add(TTree * tree, const char * branch_name,  TEventList * list, double interpolation_dt) ; 

        /*** Reset to empty state **/ 
        void reset(); 

        /*** Get the averaged waveform */ 
        TGraph * getAverage() { computeAverage(); return avg; }

        /*** Get number in average  **/ 
        int nAveraged() const { return n; }
        

      private:
        void computeAverage(); 
        TGraph * avg; 
        TGraph * sum; 
        double t0,t1; 
        double max_shift; 
        bool dirty; 
        int * norm; 
        int n;


    }; 
}


#endif 
