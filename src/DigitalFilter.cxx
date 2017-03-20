#include "DigitalFilter.h" 
#include <complex>
#include "TCanvas.h" 
#include <assert.h>
#include "TGraph.h" 
#include "TH2.h" 
#include "TAxis.h"
#include "TStyle.h" 
#include "TMath.h"
#include "FFTtools.h" 
#include "Math/Polynomial.h" 
#include "FFTWindow.h" 
#include "TEllipse.h" 

#include <strings.h>
#include "TMatrixD.h"
#include "TDecompLU.h"



void FFTtools::DigitalFilter::response(size_t n, TGraph ** amplitude_response, TGraph ** phase_response, TGraph ** group_delay) const 
{
  TGraph * gamp = 0; 
  TGraph * gph = 0; 
  TGraph * ggroup = 0; 

  if (amplitude_response)
  {
    if (*amplitude_response) 
    {
      gamp = *amplitude_response; 
      gamp->Set(n); 
   }
   else
    {
      gamp = new TGraph(n); 
      *amplitude_response = gamp; 
    }
  }

  if (phase_response)
  {
    if (*phase_response) 
    {
      gph = *phase_response; 
      gph->Set(n); 
    }
    else
    {
      gph = new TGraph(n); 
      *phase_response = gph; 
    }
     
  }

  if (group_delay) 

  {
    if (*group_delay)
    {
      ggroup = *group_delay; 
      ggroup->Set(n);
    }
    else
    {
      ggroup = new TGraph(n); 
      *group_delay = ggroup; 
    }
  }

  if (!phase_response && !amplitude_response && !group_delay) return; 

  for (size_t i = 0; i < n; i++) 
  {

    double f= 1./(n-1) *i; 
    std::complex<double> z = exp(std::complex<double>(0, TMath::Pi() * f)); 
    std::complex<double> resp = transfer(z); 
    if (gamp)
    {
      gamp->SetPoint(i, f, abs(resp) == 0 ? -1000 : 10*log(abs(resp))); 
    }
    double angle = 0; 

    if (gph || ggroup) 
    {
      angle = 180. * (arg(resp) / TMath::Pi()); 

      //round to 6th decimal place
      angle = rint(angle * 1.e6)/1.e6; 
      if (angle < 0) angle += 360;
    }

    if (gph)
    {
      gph->SetPoint(i, f, angle); 
    }
    if (ggroup) 
    {
      ggroup->SetPoint(i,f,angle) ; 
    }
  }
  if (gamp)
  {
    gamp->SetTitle("Amplitude Response"); 
    gamp->GetXaxis()->SetTitle("Normalized Frequency (f/f_{nyq})"); 
    gamp->GetYaxis()->SetTitle("Magnitude (dB)"); 
    gamp->SetMinimum(-80); 
    gamp->SetMaximum(10); 
  }
  if (gph)
  {
    gph->RemovePoint(n-1); 
    gph->RemovePoint(0);
    gph->SetTitle("Phase Response"); 
    gph->GetXaxis()->SetTitle("Normalized Frequency (f/f_{nyq})"); 
    gph->GetYaxis()->SetTitle("Phase (deg))");
  }

  if (ggroup) 
  {
    unwrap(ggroup->GetN(), ggroup->GetY(), 360); 
    double last = 2*ggroup->GetY()[0];  //TODO do I want this factor of 2? 
    ggroup->GetY()[0] = 0; 
    for (int i = 1; i < ggroup->GetN(); i++) 
    {
     double current = 2*ggroup->GetY()[i]; 
     ggroup->GetY()[i] = last - current; 
     last = current; 
    }

    ggroup->RemovePoint(n-1); 
    ggroup->RemovePoint(0);
    ggroup->RemovePoint(0);

    ggroup->SetTitle("Group Delay"); 
    ggroup->GetXaxis()->SetTitle("Normalized Frequency (f/f_{nyq})"); 
    ggroup->GetYaxis()->SetTitle("Normalized Time (t/(T))"); 
  }
}


void FFTtools::DigitalFilter::impulse(size_t n, double * out, size_t delay) const 
{
  double in[n]; 
  memset(in,0,n*sizeof(double)); 
  in [delay] = 1; 
  filterOut(n,in,out); 
}

TGraph * FFTtools::DigitalFilter::impulseGraph(size_t n, double dt, size_t delay) const 
{
  TGraph * g = new TGraph(n); 
  g->SetTitle(TString::Format("Impulse Response (T=%g, delay=%lu)", dt, delay)); 
  impulse(n,g->GetY(), delay); 
  for (size_t i = 0; i < n; i++) { g->GetX()[i] = i * dt; }
  return g; 
}

double * FFTtools::DigitalFilter::impulse(size_t n, size_t delay) const 
{
  double * out = new double[n]; 
  impulse(n,out,delay); 
  return out; 
}

void FFTtools::DigitalFilter::filterReplace(size_t n, double * w) const 
{
  double * filtered = filter(n,w); 
  memcpy(w,filtered, n * sizeof(*w)); 
  delete [] filtered; 
}

void FFTtools::DigitalFilter::filterGraph(TGraph * g, bool filterErrors) const 
{
  filterReplace(g->GetN(),g->GetY()); 

  if (g->GetEY() && filterErrors)
  {
    double *ey2 = new double[g->GetN()]; 
    for (int i = 0; i < g->GetN(); i++) ey2[i] = g->GetEY()[i]*g->GetEY()[i]; 
    filterOut(g->GetN(), ey2,g->GetEY()); 
    delete [] ey2; 
    for (int i = 0; i < g->GetN(); i++) g->GetEY()[i] = sqrt(g->GetEY()[i]); 
  }
}

double * FFTtools::DigitalFilter::filter(size_t n, const double * w) const 
{
  double * out = new double[n]; 
  filterOut(n,w,out); 
  return out; 
}


std::complex<double> FFTtools::DigitalFilterSeries::transfer(std::complex<double> z) const 
{
  std::complex<double> answer(1,0); 
  for (size_t i = 0; i < series.size(); i++) 
  {
    answer *= series[i]->transfer(z); 
  }
  return answer; 
}


void FFTtools::DigitalFilterSeries::filterOut(size_t n, const double * w, double * out) const 
{
  if (!series.size())
  {
    memcpy(out,w, n * sizeof(double)); 
    return; 
  }
  double * y= (double*) calloc(n,sizeof(double)); 

  double * newx= series.size() > 1 ? (double*) calloc(n,sizeof(double)) : 0; 
  const double * xin = w; 

 for (size_t i = 0; i < series.size();i++)
 {

    if (i > 0) 
    {
      double * temp = newx;  //store y as newx , reuse space of newx; 
      newx = y; 
      y = temp; 
      if ( i  > 1) 
      {
        //need to rezero
        memset(y, 0,n*sizeof(double)); 
      }
      xin = newx; 
    }
    series[i]->filterOut(n,xin,y); 
  }

  free(newx); 
  memcpy(out,y, sizeof(double) * n); 
  free(y); 
  return; 
}





//almost equivalent to matlab poly. (except indices are reversed so that the index matches the order)
static void poly(unsigned int n, const std::complex<double> * zeroes, std::complex<double>  * coeffs)
{


  //Start out with (x - zero_0)
  coeffs[0] = -zeroes[0]; 
  coeffs[1] = std::complex<double>(1,0); 


  for (unsigned int i = 1; i < n; i++)
  {

    //Multiply  by each (x - zero_i)
    for (unsigned int j = i+1;;  j--)
    {
      coeffs[j] = (j > 0 ? coeffs[j-1] : 0) - (j < i+1 ? coeffs[j] * zeroes[i] : 0); 
      if (j == 0) break; 
    }
  }
}

// see helpful http://octave.sourceforge.net/signal/function/sftrans.html

void FFTtools::TransformedZPKFilter::transform(FilterTopology type , double w, double dw)
{

  double W=tan(TMath::Pi()*w/2);   // prewarp frequencies 
  double Wh = tan(TMath::Pi()*(w+dw)/2);  
  double Wl = tan(TMath::Pi()*(w-dw)/2);  
  double dW = (Wh - Wl)/2 ; 

  std::vector<std::complex<double> > new_zeroes; 
  std::vector<std::complex<double> > new_poles; 
  std::complex<double> new_gain = gain; 

  if (type == LOWPASS) 
  {
    for (size_t i = 0; i < zeroes.size(); i++) 
    {
      new_zeroes.push_back(zeroes[i] * W); 
      new_gain /=W; 
    }

    for (size_t i = 0; i < poles.size(); i++) 
    {
      new_poles.push_back(poles[i] * W); 
      new_gain *=W; 
    }
  }
  else if (type == HIGHPASS)
  {
    for (size_t i = 0; i < zeroes.size(); i++) 
    {
      new_zeroes.push_back(W/zeroes[i]); 
      new_gain *=-zeroes[i]; 
    }

    //extend zeroes to poles, if necessary
    for (size_t i = zeroes.size(); i < poles.size(); i++) 
    {
      new_zeroes.push_back(0); 
    }
    

    for (size_t i = 0; i < poles.size(); i++) 
    {
      new_poles.push_back(W/poles[i]); 
      new_gain *= -1./poles[i]; 
    }

    //extend poles to zeroes, if necessary
    for (size_t i = poles.size(); i < zeroes.size(); i++) 
    {
      new_poles.push_back(0); 
    }
 
  }
  else
  {
    assert(dW); 
    if (type == BANDPASS)
    {
      for (size_t i = 0; i < zeroes.size(); i++) 
      {
        std::complex<double> b  = zeroes[i] * dW; 
        std::complex<double> x = sqrt(b*b - Wh*Wl); 
        new_zeroes.push_back(b + x); 
        new_zeroes.push_back(b - x); 
        new_gain *= 1./(2*dW); 
      }
      for (size_t i = zeroes.size(); i < poles.size(); i++) 
      {
        new_zeroes.push_back(0); 
      }

      for (size_t i = 0; i < poles.size(); i++) 
      {
        std::complex<double> b  = poles[i] * dW; 
        std::complex<double> x = sqrt(b*b - Wh*Wl); 
        new_poles.push_back(b + x); 
        new_poles.push_back(b - x); 
        new_gain *= 2*dW; 
      }

      for (size_t i = poles.size(); i < zeroes.size(); i++) 
      {
        new_poles.push_back(0); 
      }
    }
    else if (type == NOTCH)
    {
      std::complex<double> extra = sqrt(std::complex<double>(-Wh*Wl)); 

      for (size_t i = 0; i < zeroes.size(); i++) 
      {
        std::complex<double> b  =dW / zeroes[i]; 
        std::complex<double> x = sqrt(b*b - Wh*Wl); 
        new_zeroes.push_back(b + x); 
        new_zeroes.push_back(b - x); 
        new_gain *= -zeroes[i]; 
      }
      for (size_t i = zeroes.size(); i < poles.size(); i++) 
      {
        new_zeroes.push_back(extra); 
        new_zeroes.push_back(-extra); 
      }

      for (size_t i = 0; i < poles.size(); i++) 
      {
        std::complex<double> b  = dW / poles[i]; 
        std::complex<double> x = sqrt(b*b - Wh*Wl); 
        new_poles.push_back(b + x); 
        new_poles.push_back(b - x); 
        new_gain *= -1./poles[i]; 
      }

      for (size_t i = poles.size(); i < zeroes.size(); i++) 
      {
        new_poles.push_back(extra); 
        new_poles.push_back(-extra); 
      }
    }
  }


  gain = std::real(new_gain); 
  zeroes = new_zeroes; 
  poles = new_poles; 


}


void FFTtools::TransformedZPKFilter::bilinearTransform()
{
  
  size_t n = std::max(poles.size(), zeroes.size()) ; 


  digi_poles.clear(); digi_poles.insert(digi_poles.begin(), n,0); 
  digi_zeroes.clear(); digi_zeroes.insert(digi_zeroes.begin(), n,0); 


  digi_gain = std::complex<double>(gain,0);  

  std::complex<double> one(1,0); 
  for (unsigned int i = 0; i < n; i++)
  {

    if (i < zeroes.size())
    {
      digi_gain *=  (one -zeroes[i]); 
    }
   
    if (i < poles.size())
    {
      digi_gain /= (one - poles[i]); 
    }

    digi_poles[i] = i < poles.size() ? (one + poles[i]) / (one - poles[i]) : std::complex<double>(-1,0); 
    digi_zeroes[i] = i < zeroes.size() ? (one + zeroes[i]) / (one - zeroes[i]) : std::complex<double>(-1,0); 
  }

  digi_gain = std::real(digi_gain); 

  //Now, get coefficients from poles and zeroes

  computeCoeffsFromDigiPoles(); 
}


void FFTtools::FIRFilter::filterOut(size_t n, const double* x, double * out) const
{
  directConvolve(n,x,coeffs.size(), &coeffs[0], out, delay, extend ? REPEAT_OUTSIDE : ZEROES_OUTSIDE); 
}


static void computeSavitzkyGolayCoefficients(double *x, int m, int nl, int nr, int deriv = 0) 
{
  assert(nl >= 0 && nr >= 0 && m >= 0 && m > deriv); 
  int size = nl+nr+1; 

  //just use ROOT TMatrix here... performance is not crucial 
  TMatrixD M(size,m+1); 
  for (int i = 0; i < size; i++) 
  {
    for (int j = 0; j <= m; j++)
    {
      M(i,j) = TMath::Power(i-nl, j); 
    }
  }
  TMatrixD Mcopy(M); 
  TDecompLU lu(Mcopy.T() * M); 
  TVectorD y(m+1); 
  y(deriv) = 1; 
  lu.Solve(y); 


  double factorial = TMath::Factorial(deriv); 
  for (int i = 0; i < size; i++) 
  {
    double dotproduct = y(0); 
    for (int j = 1; j <=m; j++) 
    {
      dotproduct += y(j) * M(i,j);
    }
    x[i] = dotproduct*factorial; 
  }
}



FFTtools::SavitzkyGolayFilter::SavitzkyGolayFilter(int polynomial_order, int wleft, int wright, int deriv) 
  : FIRFilter(wright < 0 ? 2*wleft + 1 : wleft +wright +1, wright < 0 ? 0 : wright- wleft,true)
{
  if (wright < 0) wright = wleft; 
  computeSavitzkyGolayCoefficients(&coeffs[0], polynomial_order, wleft, wright, deriv); 
}


std::complex<double> FFTtools::FIRFilter::transfer(std::complex<double> z) const
{
  std::complex<double> ans = 0; 
  for( size_t i = 0; i < coeffs.size(); i++) 
  {
    int exp = int(coeffs.size()/2) - int(i) - delay; 
    ans += pow(z,exp) * coeffs[i]; 
  }

  return ans; 

}

FFTtools::SincFilter::SincFilter(double w, int max_lobes, const FFTWindowType * win, int delay, bool extend) 
  : FIRFilter(2*max_lobes/w + 1, delay, extend)

{
  int N = 2 * max_lobes / w + 1; 
  for (int i = 0; i < N; i++)  
  {
   coeffs[i] =  w * sinc (w * (i - N/2)); 
  }

  if (win) win->apply(coeffs.size(),&coeffs[0]);
}

static void doIIRFilter(int n, const double * x, double * y, int na, const double * A, int nb, const double * B)
{
  int nc = TMath::Max(na,nb); 
  for (int j = 0; j < n; j++) 
  {
      double a0 =A[0];
      y[j] = 0; 
      for (int k = 0; k < nc; k++) 
      {
        if (j - k < 0) break; 

        if (k < nb)
        {
          y[j] += x[j-k] * B[k] ; 
        }

        if (k > 0 && k < na) 
        {
          y[j] -= y[j-k] *A[k]; 
        }
      }
      y[j] /= a0;
  }
}


void FFTtools::IIRFilter::filterOut(size_t n, const double* x, double * out) const
{
  doIIRFilter(n, x, out, acoeffs.size(), &acoeffs[0],bcoeffs.size(), &bcoeffs[0]); 
}

std::complex<double> FFTtools::IIRFilter::transfer(std::complex<double> z) const
{
  std::complex<double> num = bcoeffs[0]; 
  for( size_t i = 1; i < bcoeffs.size(); i++) 
  {
    num += std::pow(z,-double(i)) * bcoeffs[i]; 
  }

  std::complex<double> denom = acoeffs[0]; 
  for( size_t i = 1; i < acoeffs.size(); i++) 
  {
    denom += std::pow(z,-double(i)) * acoeffs[i]; 
  }

  std::complex<double> answer = num/denom; 

  return answer; 
}


FFTtools::RCFilter::RCFilter(FilterTopology type, double w, double dw)
{
  order = 1; 
  poles.push_back(-1); 
  gain = 1; 
  transform(type,w,dw); 
  bilinearTransform(); 
}

FFTtools::ButterworthFilter::ButterworthFilter(FilterTopology type, size_t order,  double w, double dw )
{
  this->order = order; 
  poles.reserve(order);

  std::complex<double> pole_product(1,0); 

  //Compute Butterworth poles for given order
  for (unsigned int i = 0; i < order; i++)
  {
    poles.push_back(std::exp(std::complex<double>(0, TMath::Pi()* (2*(i+1)+order - 1)/(2*order)))); 
    pole_product*= -poles[i]; 
  }

  gain = 1./std::real(pole_product); 

  transform(type,w,dw); 
  bilinearTransform(); 
}


FFTtools::ChebyshevIFilter::ChebyshevIFilter(FilterTopology type, size_t order,  double ripple, double w, double dw )
{
  this->order = order; 
  poles.reserve(order);

  std::complex<double> pole_product(1,0); 
  double eps = sqrt(TMath::Power(10,ripple/10) -1); 
  double v0 = TMath::ASinH(1/eps)/order; 
  double sv0 = TMath::SinH(v0); 
  double cv0 = TMath::CosH(v0); 


  //Compute Cheby1 poles for given order
  for (unsigned int i = 0; i < order; i++)
  {

    double theta = TMath::Pi() / 2 * (2*(i+1)-1) / order; 
    poles.push_back(std::complex<double>(-sv0 * sin(theta),cv0*cos(theta))); 
    pole_product*= -poles[i]; 
  }

  gain = std::real(pole_product); 

  if (order % 2 == 0) gain/=TMath::Power(10,ripple/20); //adjust for ripple  start

  transform(type,w,dw); 
  bilinearTransform(); 
}

FFTtools::GaussianFilter::GaussianFilter(double sigma, double nsigma) 
  : FIRFilter (2* ceil(nsigma * sigma) + 1)
{
  int w = ceil(nsigma * sigma); 
  double sum = 0; 
  for (int i = 0; i < 2 * w+ 1; i++) 
  {
    double x = (i - w); 
    coeffs[i] = TMath::Gaus(x,0,sigma); 
    sum += coeffs[i]; 
  }

  for (int i = 0; i < 2 * w+ 1; i++) 
  {
    coeffs[i] /= sum; 
  }
 
}

FFTtools::BoxFilter::BoxFilter(int width)
  : FIRFilter (width) 
{
  for (int i = 0; i < width; i++)
  {
    coeffs[i] = 1./width; 
  }
}

FFTtools::DifferenceFilter::DifferenceFilter(int order)
  : FIRFilter (order +1) 
{
  assert(order > 0); 
  double sum = 0; 
  for (int i = 0; i <= order; i++)
  {
    coeffs[i] = TMath::Binomial(order,i)  * (i %2 == 0 ? 1 : -1); 
    sum += fabs(coeffs[i]); 
  }


  for (int i = 0; i <= order; i++)
  {
    coeffs[i] /= sum; 
  }
}

void FFTtools::IIRFilter::computePolesAndZeroes() 
{
  digi_gain = bcoeffs[0]; 

  ROOT::Math::Polynomial pb(bcoeffs.size()-1); 
  std::vector<double> scaled_bcoeffs(bcoeffs.size()); 

  for (unsigned i = 0; i < bcoeffs.size(); i++) 
  {
    scaled_bcoeffs[i] = bcoeffs[i] / bcoeffs[0]; 
  }

  pb.SetParameters(&scaled_bcoeffs[0]); 
  digi_zeroes = pb.FindRoots(); 
}

FFTtools::IIRFilter::IIRFilter(const FIRFilter & other) 
  : order(other.getOrder() + other.getDelay()) , bcoeffs(order), gzero(0), gpole(0) 
{
  memcpy(&bcoeffs[other.getDelay()], other.getCoeffs(), other.getOrder()*sizeof(double)); 
  acoeffs.push_back(1); 
  computePolesAndZeroes(); 
}



void FFTtools::IIRFilter::computeCoeffsFromDigiPoles()
{

  unsigned nzeroes = nDigiZeroes(); 
  unsigned npoles = nDigiPoles(); 

  bcoeffs.clear(); bcoeffs.insert(bcoeffs.begin(), nzeroes+1,0); 
  acoeffs.clear(); acoeffs.insert(acoeffs.begin(), npoles+1,0); 

  std::vector< std::complex<double> >bpoly(nzeroes+1);
  poly(nzeroes, &digi_zeroes[0],&bpoly[0]); 

  std::vector<std::complex<double> >apoly(npoles+1);
  poly(npoles, &digi_poles[0],&apoly[0]); 

  for (unsigned int i = 0; i < nzeroes+ 1; i++)
  {
       bcoeffs[i] = (double) std::real(digi_gain * bpoly[nzeroes - i]);             
  }

  for (unsigned int i = 0; i < npoles+ 1; i++)
  {
       acoeffs[i] = (double) std::real(apoly[npoles - i]);             
  }
}

static TEllipse ell(0,0,1,1); 

void FFTtools::IIRFilter::Draw(const char * opt, int color) const
{
  bool same = strcasestr(opt,"Same"); 
//  bool pol = !same && strcasestr(opt,"pol"); 

  
  if (!gpole) 
  {
    gpole = new TGraph(nDigiPoles()); 
    for (size_t i = 0; i < nDigiPoles(); i++) 
    {
      gpole->SetPoint(i, digi_poles[i].real(), digi_poles[i].imag()); 
    }
    gpole->SetMarkerStyle(5); 

  }

  if (!gzero) 
  {
    gzero = new TGraph(nDigiZeroes()); 
    for (size_t i = 0; i < nDigiZeroes(); i++) 
    {
      gzero->SetPoint(i, digi_zeroes[i].real(), digi_zeroes[i].imag()); 
    }
    gzero->SetMarkerStyle(4); 
  }



  if (!same) 
  {
    /* figure out  maximum abs */ 

    double max = 0; 

    for (size_t i = 0; i < nDigiPoles(); i++) 
    {
      if (std::abs(digi_poles[i]) > max)
        max = std::abs(digi_poles[i]); 
    }

    for (size_t i = 0; i < nDigiZeroes(); i++) 
    {
      if (std::abs(digi_zeroes[i]) > max)
        max = std::abs(digi_zeroes[i]); 
    }

    int optstat = gStyle->GetOptStat(); 
    int optfit = gStyle->GetOptFit(); 
    TH2I ax("polaxis","Pole-Zero Plot", 10,-max,max,10,-max,max); 
    gStyle->SetOptStat(0); 
    ax.DrawCopy("a"); 
    gStyle->SetOptStat(optstat); 
    gStyle->SetOptFit(optfit); 
    ell.Draw(); 
  }


  gpole->SetMarkerColor(color); 
  gzero->SetMarkerColor(color); 

  gpole->Draw("msame"); 
  gzero->Draw("msame"); 
}




static void reflectRoots(int N, const std::complex<double> *  roots, std::vector<std::complex<double> > & output, double eps) 
{


  for (int i = 0; i < N; i++) 
  {
    double M = abs(roots[i]); 
    if ( M < 1) // woohoo 
    {
      output[i] = roots[i]; 
    }
    else if (M > 1)  //reflect it 
    {
      output[i] = 1./roots[i]; 
    }
    else //attenuate it a bit
    {
      output[i] = roots[i] - std::complex<double>(eps,eps); 
    }
  }
}




FFTtools::MinimumPhaseFilter::MinimumPhaseFilter(const IIRFilter & other, double eps) 
{

  
  order = other.getOrder(); 
  digi_gain = other.getDigiGain();
  reflectRoots(other.nDigiPoles(), other.getDigiPoles(), digi_poles, eps); 
  reflectRoots(other.nDigiZeroes(), other.getDigiZeroes(), digi_zeroes, eps); 

  computeCoeffsFromDigiPoles(); 
}










TString FFTtools::IIRFilter::asString() const 
{

  /*This is woefully inefficient, but hopefully nobody uses this in a loop */ 


  TString top =  "      "; 
  TString mid =  "H(z)= "; 
  TString bot =  "      "; 

  top += TString::Format("%.4g ",bcoeffs[0]); 

  for (unsigned i = 1; i < bcoeffs.size(); i++) 
  {
    top += TString::Format(" + %.4g z^-%u", bcoeffs[i],i); 
  }


  bot += TString::Format("%.4g ",acoeffs[0]); 
  for (unsigned i = 1; i < acoeffs.size(); i++) 
  {
    bot += TString::Format(" + %.4g z^-%u", acoeffs[i],i); 
  }

  if (top.Length() < bot.Length())
  {
    int diffmid = bot.Length() - mid.Length(); 
    for (int i = 0; i < diffmid; i++) mid+= "-"; 
    int diff = bot.Length()-top.Length(); 
    for (int i = 0; i < (diff)/2; i++) top = " " + top; 
  }
  else
  {

    int diffmid = top.Length() - mid.Length(); 
    for (int i = 0; i < diffmid; i++) mid+= "-"; 
    int diff = top.Length()-bot.Length(); 
    for (int i = 0; i < (diff)/2; i++) bot = " " + bot; 
  }

  return "IIR Transfer Function: \n" + top + "\n" + mid + "\n" + bot + "\n"; 
}


FFTtools::ThiranFilter::ThiranFilter(double delay, int order)
{
  //implementation based on https://ccrma.stanford.edu/~jos/pasp/Thiran_Allpass_Interpolators.html

  if (delay < 0) 
  {

    fprintf(stderr,"ThiranFilter requires non-negative delay\n"); 
    order = 0; 
  }

  if (order < 1) 
  {
    order = ceil(delay); 
  }

  int N = order; //to make code more readable 

  acoeffs.push_back(1); 
  for (int k = 1; k <= N; k++) 
  {
    double a_k = ( (k % 2)  ? -1 : 1) * TMath::Binomial(N,k); 

    for (int n = 0; n <= N; n++)
    {
      a_k *= (delay - N + n); 
      a_k /= (delay - N + k + n); 
    }

    acoeffs.push_back(a_k); 
  }

  //the bcoeffs will be the acoeffs backwards

  for (int i = N; i >=0; i--) 
  {
    bcoeffs.push_back(acoeffs[i]); 
  }

}



TPad* FFTtools::DigitalFilter::drawResponse( TPad * c, int n, int delay)  const
{
  if (!c) c = new TCanvas("filterResponse","Filter Response"); 
  c->Divide(2,2);; 

  TGraph * amp = 0; 
  TGraph * ph = 0; 
  TGraph * grp = 0; 
  response(n,&amp,&ph,&grp); 

  c->cd(1); 
  amp->Draw(); 
  c->cd(2); 
  ph->Draw(); 
  c->cd(3); 
  grp->Draw(); 
  c->cd(4); 
  impulseGraph(n,1,delay)->Draw("alp"); 
 

  return c; 

}
