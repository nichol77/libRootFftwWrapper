#ifndef _FFTTOOLS_CWT_HH
#define _FFTTOOLS_CWT_HH


/* Continuous wavelet transform  implementation 
 *
 *  Cosmin Deaconu <cozzyd@kicp.uchicago.edu> 
 * */ 



#include <cstdlib>
#include <vector>

class TGraph; 
class TH2; 


namespace FFTtools
{


    class CWT
    {
          public: 

            /* Perhaps templates would be more performant here to avoid the virtual table overhead */ 
            class MotherWavelet
            {
              public: 
                virtual void fill(int M, double width, double *wavelet) const; 
                virtual double eval(double x, double w) const = 0; 
                virtual ~MotherWavelet() {;}
            }; 

            class Ricker : public MotherWavelet 
            {
              public:
                virtual void fill(int M, double width, double * wavelet) const ; 
                virtual double eval(double x, double w) const; 
                virtual ~Ricker() {; }
            };

            class Ridger : public MotherWavelet 
            {
              public:
                void fill(int M, double width, double * wavelet) const; 
                virtual double eval(double x, double w) const; 
                virtual ~Ridger() {; }
            };

            class Morlet : public MotherWavelet 
            {
              public:
                Morlet(double freq = 5) { f = freq; }
                void fill(int M, double width, double * wavelet) const; 
                virtual double eval(double x, double w) const; 
                virtual ~Morlet() {; }
              private:
                double f; 
            };



       
           
        CWT(const MotherWavelet *w, double scaled_width =5, int nscales =10 , double start_scale = 1, double scale_step = 1); 
        virtual ~CWT() { cleanup(); } 

        /* Compute CWT; assume that TGraph has uniform spacing (has been interpolated) */ 
        void compute(const TGraph * win); 

        void setScales(int n, const double * scales); 
        void setScales(int n, double start_scale, double scale_step); 
        void setWidth(double w) { sw = w; }

        void setMotherWavelet(const MotherWavelet * w) {mw = w;}

        size_t nScales() const { return scales.size(); }
        const double *getScales() const { return &scales[0]; }
        double getScale(int scale) { return scales[scale]; }

        double eval(double t, double min_abs_val= 0, int min_scale = 0, int max_scale = 0) const; 
        void eval(int N,  const double *t, double * out, double min_abs_val= 0, int min_scale = 0, int max_scale = 0) const; 



        size_t getSize() const { return size; }
        const double * getX() const { return xvec; } 
        const double * getY(int scale) const { return yvecs[scale]; } 
        TGraph * getGraph(int scale) const; 
        TH2 * getHist(const char * name = "cwt") const ; 

        

      private: 

        const MotherWavelet * mw; 
        double sw; 
        std::vector<double> scales; 

        size_t size; 
        double * xvec; 
        std::vector<double *> yvecs; 
        void cleanup(); 
        

    };

}





#endif

