#include "TGraph.h" 
#include "TGraphErrors.h" 
#include "RFInterpolate.h" 
#include "TH2.h"
#include <iostream>
#include "FFTtools.h" 
#include <assert.h>

#ifdef USE_EIGEN
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseQR>

#define MAT Eigen::MatrixXd 
#define SPMAT Eigen::SparseMatrix<double> 
#define VEC Eigen::VectorXd 

// bdcSvd only in Eigen 3.3+, much faster than jacobiSvd 
//#define LINSOLVE(A, B)  A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(B)
//
//Use QR decomposition 
#define LINSOLVE(A, B)  A.householderQr().solve(B)

//Use Normal Equations 
//#define LINSOLVE(A,B) (A.transpose() *A).ldlt().solve(A.transpose() *B) 

#define RESIZE_MAT(A, m,n) A.conservativeResize(m,n) 
#define RESIZE_VEC(A, n) A.conservativeResize(n) 
#define NROWS(A) A.rows() 
#define NCOLS(A) A.cols() 

#elif USE_ARMADILLO
#include <armadillo> 
#define MAT arma::mat
#define SPMAT arma::SpMat<double>
#define VEC arma::vec

#define LINSOLVE(A, B)  arma::solve(A,B) 
#define RESIZE_MAT(A, m,n) A.set_size(m,n) 
#define RESIZE_VEC(A, n) A.set_size(n)
#define NROWS(A) A.n_rows 
#define NCOLS(A) A.n_cols 

#else //use slow ROOT methods
#include "TMatrixD.h"
#include "TVectorD.h"
//#include "TDecompInvert.h"
//#include "TDecompSVD.h"
#include "TDecompQRH.h"
#include "TDecompSparse.h"

#define MAT TMatrixD
#define VEC TVectorD
#define SPMAT TMatrixDSparse 

#define RESIZE_MAT(A, m,n) A.ResizeTo(m,n) 
#define RESIZE_VEC(A, n) A.ResizeTo(n) 
#define NROWS(A) A.GetNrows()  
#define NCOLS(A) A.GetNcols() 

static bool ignore_me_at_your_peril;  
#define LINSOLVE(A,B)  TDecompQRH(A).Solve(B,ignore_me_at_your_peril) 
#endif 



#include "TMath.h" 
//#include "TStopwatch.h" 

/*
static std::complex<double> sinc(std::complex<double> x)
{
  if (x == 0) 
  {
    return 1; 
  }

  return sin(TMath::Pi() * x) / (TMath::Pi() * x); 
}
*/

// d/(dx) (sin pi x) / (pi x) =  cos(pi * x) / x - sin(pi * x) / (pi * x^2) 
static double sincprime (double x, double eps = 0) 
{

  if (fabs(x) <= eps)
  {
    return 0; 
  }
  else return cos(TMath::Pi() * x)/x - FFTtools::sinc(x)/x; 

}

static void infer_vals(double & dt, int n, int & nout, const double * xj)
{
  bool given_dt = dt >=0; 
  if (dt <= 0) 
  {
    dt = (xj[n-1] - xj[0])/(n-1); 
  }

  if (nout <= 0)
  {
    nout  = given_dt ?  ceil((xj[n-1]-xj[0])/dt) : n; 
  }
}


TGraph * FFTtools::getInterpolatedGraphInvertSplit(const TGraph * g, double dt, int size, int overlap, int nout) 
{
//  TStopwatch w; 

  const double * xj = g->GetX(); 
  const double * yj = g->GetY(); 
  int n=  g->GetN();
  
  infer_vals(dt,n,nout,xj); 
  double t0 = xj[0]; 
  int nsegs = ceil(double(nout) / size); 


  MAT A(size+overlap,size+overlap);
  VEC B(size+overlap);

  double x[nout]; 
  double y[nout]; 

  for (int i = 0; i < nout; i++) 
  {
    x[i] = t0 + i *dt; 
    y[i] = 0; 
  }


  int start = 0; 


  for (int seg = 0; seg < nsegs; seg++) 
  {
    int end = start + size > nout ? nout: start+size; 
    double start_out = start *dt + t0;  
    double end_out = end *dt + t0;  
    int in_start = TMath::BinarySearch(n, xj, start_out) - overlap; 
    int in_end = 1+TMath::BinarySearch(n, xj, end_out) + overlap; 
    if (in_start < 0) in_start = 0; 
    if (in_end > n) in_end = n; 

    int size_in = in_end - in_start; 
    int size_out = end - start; 

    RESIZE_MAT(A,size_in,size_out);
    RESIZE_VEC(B,size_in);

    for (int i = 0; i < size_in; i++) 
    {
      B(i) = yj[in_start + i]; 
      for (int j = 0; j < size_out; j++) 
      {
        A(i,j) = sinc((x[j+start] - xj[i+in_start])/dt); 
      }
    }

    VEC soln = LINSOLVE(A,B);

    for (int i = 0; i < size_out; i++) 
    {
      y[i+start] = soln(i); 
    }
    start += size; 
  }


  TGraph * gout =  new TGraph(nout,x,y); 
//  w.Print(); 
  return gout; 
}

TGraph * FFTtools::getInterpolatedGraphInvertLapped(const TGraph * g, double dt, int size, int nout) 
{
//  TStopwatch w; 

  const double * xj = g->GetX(); 
  const double * yj = g->GetY(); 
  int n =  g->GetN();

  infer_vals(dt,n,nout,xj); 

  assert(size % 2 == 0); 

  double t0 = xj[0]; 

  int nsegs = ceil(2.*nout / size -1); 



  double x[nout]; 
  double y[nout]; 

  for (int i = 0; i < nout; i++) 
  {
    x[i] = t0 + i *dt; 
    y[i] = 0; 
  }


  int start = 0; 

  MAT A(size,size); 
  VEC B(size); 

  for (int seg = 0; seg < nsegs; seg++) 
  {

    int end = start + size > nout ? nout : start+size; 

    double start_out = start *dt + t0;  
    double end_out = end *dt + t0;  
    int in_start = TMath::BinarySearch(n, xj, start_out); 
    int in_end = 1+TMath::BinarySearch(n, xj, end_out); 
    if (in_end > n) in_end = n; 

    int size_in = in_end - in_start; 
    int size_out = end - start; 


    RESIZE_MAT(A,size_in,size_out);
    RESIZE_VEC(B,size_in);

    for (int i = 0; i < size_in; i++) 
    {
      B(i) = yj[in_start + i]; 


      for (int j = 0; j < size_out; j++) 
      {
         A(i,j) = sinc((x[j+start] - xj[i+in_start])/dt); 
      }
    }

    VEC soln = LINSOLVE(A,B);

    for (int i = 0; i < size_out; i++) 
    {
      if (i + start < size/2 || i + start >= nout-size/2 ) 
      {
        y[i+start] = soln(i); 
      }
      else 
      {
        double weight = (i <= size/2) ? 2*double(i)/size  : 2-2*double(i)/size; 
        y[i+start] += soln(i)*weight;  
      }
    }
    start += size/2; 
  }


  TGraph * gout =  new TGraph(nout,x,y); 
//  w.Print(); 
  return gout; 
}

#ifdef USE_EIGEN
//static Eigen::SparseQR<Eigen::SparseMatrix<double> , Eigen::COLAMDOrdering<int> > solver; 
static Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver; 
//static Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > solver; 
#endif 

TGraphErrors * FFTtools::getInterpolatedGraphSparseInvert(const TGraph * g, double dt, int nout, double max_dist, 
                                                        double eps, double weight_exp, double lambda, int lambda_order, double error_scale, 
                                                         TH2 * hA) 
{

#if !defined(USE_EIGEN) && !defined(USE_ARMADILLO) 
  static bool already_scolded = false; 
  if (!already_scolded) printf("WARNING: Strongly recommend compiling with eigen3 (or even armadillo) support when calling getInterpolatedGraphSparseInvert. See Makefile.config in libRootFftwWrapper. \n"); 
  already_scolded = true; 
#endif 

  const double * xj = g->GetX(); 
  const double * yj = g->GetY(); 
  int n =  g->GetN();
  infer_vals(dt,n,nout,xj); 
  double t0 = xj[0]; 

  double x[nout]; 
  double y[nout]; 
  double ey[nout]; 
  int ny[nout]; 


  for (int i = 0; i < nout; i++) 
  {
    x[i] = t0 + i *dt; 
    ey[i] =0; 
    ny[i] =0; 
  }



  SPMAT A(n,nout); 
  VEC B(n); 

#ifdef USE_EIGEN 
  std::vector<Eigen::Triplet<double> > triplets; 
  triplets.reserve(2*max_dist / dt*n); 
#endif 

  if (hA) 
  {
    hA->Reset(); 
    hA->SetBins(nout,0,nout,n,0,n); 
  }

  double avg_dt = xj[n-1]/(n-1); 
  for (int i = 0; i < n; i++) 
  {
    double weight = g->GetEY() ? 1./(g->GetEY()[i] * g->GetEY()[i]) :  1.; 

    if (weight_exp > 0)
    {
      double this_dt = i ==0 ? dt : xj[i] - xj[i-1]; 
      weight *= TMath::Power(avg_dt/this_dt,weight_exp); 
    }

    B(i) = weight * yj[i]; 
    for (int j = 0; j < nout; j++) 
    {
      if (max_dist <=0 || fabs(x[j] - xj[i]) <=max_dist *dt)
      {

         double dx = (x[j]-xj[i])/dt; 
         double val = sinc(dx); 
         ey[j] +=weight/(1+dx*dx); 
         ny[j] ++; 
         val *= weight; 
//         printf("%f %f %d %d %f\n", val, weight, i,j, xj[i]-x[j]); 
         if (fabs(val) <= eps) 
         {
           continue; 
         }
#ifdef USE_EIGEN
         triplets.push_back(Eigen::Triplet<double>(i, j, val)); 
#else
         A(i,j) = val; 
#endif
         if (hA)
         {
           hA->SetBinContent(1+j,n-i,val); 
         }
      }
    }
  }




#ifdef USE_EIGEN
  A.setFromTriplets(triplets.begin(), triplets.end()); 


//  std::cout << A << std::endl; 
//  solver.compute(A); 
//  VEC soln = solver.solve(B); 
  SPMAT id(nout,nout); 
  if (lambda_order == 0)
  {
    id.setIdentity(); 
  }
  else
  {
    SPMAT tmp(nout-lambda_order,nout); 
    triplets.clear(); 
    for (int i = 0; i < nout-lambda_order; i++) 
    {
      for (int j = 0; j <= lambda_order; j++) 
      {
        triplets.push_back(Eigen::Triplet<double>(i,i+j, TMath::Binomial(lambda_order,j) * (j %2 == 0 ? -1 : 1)));
      }
    }

    tmp.setFromTriplets(triplets.begin(), triplets.end()); 
    id = (tmp.transpose() * tmp); 
  }

  SPMAT AtA = A.transpose()*A; 

  double trace = 0; 
  for (int i = 0; i < nout; i++) 
  {
    trace += AtA.coeff(i,i); 
  }

//  printf("%f %f\n",trace, lambda*trace/nout); 
  solver.compute(AtA+ (lambda * trace/nout)*id); 
  VEC soln = solver.solve(A.transpose()*B); 
#elif USE_ARMADILLO
  
  VEC soln = arma::spsolve(A,B); 
#else 
  VEC soln; 
  if (n == nout) 
  {
    TDecompSparse solver(A,0); 
    soln = solver.Solve(B,ignore_me_at_your_peril); 
  }
  else
  {
    TDecompSparse solver(A.T()*A,0); 
    soln = solver.Solve(A.T() *B,ignore_me_at_your_peril); 

  }
#endif 


  for (int i = 0; i < nout; i++)
  {
    y[i] = soln(i); 
    ey[i] = ny[i] ? error_scale/sqrt(ey[i]/ny[i]) : TMath::Infinity() ; 
  }

  return new TGraphErrors(nout,x,y,0,ey); 
}


TGraph * FFTtools::getInterpolatedGraphWeightedInvert(const TGraph * g, double dt, int nout) 
{
  const double * xj = g->GetX(); 
  const double * yj = g->GetY(); 
  int n =  g->GetN();
  infer_vals(dt,n,nout,xj); 
  double t0 = xj[0]; 

  double x[nout]; 
  double y[nout]; 


  for (int i = 0; i < nout; i++) 
  {
    x[i] = t0 + i *dt; 
  }



  MAT A(n,nout); 
  VEC B(n); 

  double last_dt = dt; 
  for (int i = 0; i < n; i++) 
  {
    double this_dt = i == n-1 ? dt : xj[i+1] - xj[i]; 
    double avg_dt = (this_dt + last_dt)/2; 
    last_dt = this_dt; 
//    double weight = TMath::Power(dt/(avg_dt),2);  
    double weight = dt/avg_dt; 

    B(i) = weight*yj[i]; 
    for (int j = 0; j < nout; j++) 
    {
       A(i,j) = weight*sinc((x[j] - xj[i])/dt); 
    }
  }

  VEC soln = LINSOLVE(A,B); 


  for (int i = 0; i < nout; i++)
  {
    y[i] = soln(i); 
  }

  return  new TGraph(nout,x,y); 
}
TGraph * FFTtools::getInterpolatedGraphInvert(const TGraph * g, double dt, int nout) 
{
  const double * xj = g->GetX(); 
  const double * yj = g->GetY(); 
  int n =  g->GetN();
  infer_vals(dt,n,nout,xj); 
  double t0 = xj[0]; 

  double x[nout]; 
  double y[nout]; 


  for (int i = 0; i < nout; i++) 
  {
    x[i] = t0 + i *dt; 
  }



  MAT A(n,nout); 
  VEC B(n); 

  for (int i = 0; i < n; i++) 
  {
    B(i) = yj[i]; 
    for (int j = 0; j < nout; j++) 
    {
       A(i,j) = sinc((x[j] - xj[i])/dt); 
    }
  }

  VEC soln = LINSOLVE(A,B); 


  for (int i = 0; i < nout; i++)
  {
    y[i] = soln(i); 
  }

  return  new TGraph(nout,x,y); 
}



TGraph * FFTtools::supersample(const TGraph *g, int supersample_factor, int radius)
{

  int ncoeffs = 2*supersample_factor*radius +1; 

  //TODO: this can be done at compile time 
  double sinctable[ncoeffs]; 
  for (int i = -radius*supersample_factor; i<= radius*supersample_factor; i++) 
  {
    sinctable[radius*supersample_factor+i] = sinc(double(i)/supersample_factor); 
  }

  TGraph * newg = new TGraph( g->GetN() * supersample_factor); 
  double t0 = g->GetX()[0]; 
  double dt = (g->GetX()[1] - t0) / supersample_factor; 
  for (int i = 0; i < newg->GetN(); i++) 
  {
      double t =t0 + i*dt; 
      newg->SetPoint(i,t, (i % supersample_factor) == 0 ?  g->GetY()[i/supersample_factor] : 0);
  }

  double *convolved = directConvolve(newg->GetN(), newg->GetY(), ncoeffs,sinctable,0); 

  memcpy(newg->GetY(), convolved, newg->GetN() * sizeof(double)); 
  delete convolved; 
    

  return newg; 
}

double FFTtools::shannonWhitakerInterpolateValue(double t,  const TGraph * g, int max_lobes, const FFTWindowType * win) 
{
  const double * x = g->GetX(); 
  const double * y = g->GetY(); 
  int nx = g->GetN(); 
  double dt = x[1] - x[0]; 
  int start = 0; 
  int end = nx; 

  double center = (t-x[0])/dt; 
  int icenter = round(center); 
  if (x[icenter] == t) 
  {
    return g->GetY()[icenter]; 
  }
 
  if (max_lobes) 
  {
    start = TMath::Max(0,(int) ceil(center - max_lobes)); 
    end = TMath::Min(nx,(int) floor(center + max_lobes)); 
  }


  double val = 0; 
  for (int i = start; i < end; i++)
  {
    double arg =(t-x[i])  / dt; 
    double w = 1; 
    if (win) w = win->value(arg, max_lobes); 
    val += w*y[i]* sinc(arg); 
    
  }
  return val; 
}

double FFTtools::linearInterpolateValueAndError(double t, const TGraph* g, double * err) 
{

  const double * x = g->GetX(); 
  const double * y = g->GetY(); 
  const double * ey = g->GetEY(); 
  int nx = g->GetN(); 
  double dt = x[1] - x[0]; 

  int i = (t-x[0])/dt; 
//  printf("%d\n",i); 

  if (i < 0) 
  {
    //stupidly increase error depending on how far away it is 
    if (err) 
    {
      if (ey) 
      {
        *err = (1+20*i*i) * ey[0]; 
//        printf("%f %d %f\n", t, i, *err); 
      }
      else 
      {
        *err = 0; 
      }
    }
    return y[0]; 
  }

  else if (i >= nx-1)
  {
    if (err) 
    {
      if (ey) 
      {
        *err = (1+(i -nx)*(i-nx)) * ey[nx-1]; 
      }
      else 
      {
        *err = 0; 
      }
    }
 
    return y[nx-1]; 
  }

  double frac  = (t - x[i])/ dt; 
      
  if (err) 
  {
    if (ey)
    {
      *err = sqrt(0.5 * (frac * ey[i+1]* ey[i+1] + (1-frac) * ey[i] * ey[i])); 
    }
    else 
    {
      *err = 0; 
    }
  }

  return frac * y[i+1] + (1-frac) * y[i]; 
}

double FFTtools::shannonWhitakerInterpolateValueAndError(double t, const TGraphErrors * g, 
  double * err, int max_lobes, const FFTWindowType * win ) 
{
  const double * x = g->GetX(); 
  const double * y = g->GetY(); 
  const double * ey = g->GetEY(); 
  int nx = g->GetN(); 
  double dt = x[1] - x[0]; 

  int start = 0; 
  int end = nx; 

  double center = (t-x[0])/dt; 
  int icenter = round(center); 
  if (x[icenter] == t) 
  {
    if (err) *err = g->GetEY()[icenter]; 
    return g->GetY()[icenter]; 
  }
  if (max_lobes) 
  {
    start = TMath::Max(0,(int)ceil(center - max_lobes)); 
    end = TMath::Min(nx,(int)floor(center + max_lobes)); 
  }


 
  double val = 0; 
  double error = 0; 
//  double sum_w = 0; 
  bool nopoints = true; 
  for (int i = start; i < end; i++) 
  {
    double arg =(t-x[i])  / dt; 
    double derr = ey[i]*ey[i]; 
    double w = 1/(derr); 
    if (w == 0) continue; 
    nopoints = false; 

    double wdw = 1; 
    if (win) wdw = win->value(arg,max_lobes); 
    val += y[i]* sinc(arg)*wdw;
//    printf("%d %f %f %f\n",i, wdw,val, y[i]); 
    error += derr; 
  }

  if (nopoints)
  {
    if (err) *err =TMath::Infinity(); //infty 
    return -999; 
  }


  if (err) *err = sqrt(error); 
  return val; 
}



double FFTtools::shannonWhitakerInterpolateDerivative(double t, const TGraph * g) 
{
  const double * x = g->GetX(); 
  const double * y = g->GetY(); 
  int nx = g->GetN(); 
  double dt = x[1] - x[0]; 

  double val = 0; 
  for (int i = 0; i < nx; i++) 
  {
    double arg =(t-x[i])  / dt; 
    val += y[i]* sincprime(arg)/dt; 
  }
  return val; 
}





static const double pi = TMath::Pi(); 

TGraph * FFTtools::getInterpolatedGraphLagrange(const TGraph * g, double dt, double supersample) 
{

//  TStopwatch w; 
  const double * xj = g->GetX(); 
  int n = g->GetN(); 
  double t0 = xj[0]; 
  if (dt == 0) 
  {
    dt = (xj[n-1] - t0)/(n-1); 
  }

  double x[n]; 
  double y[n]; 


  for (int xi = 0; xi < n; xi++) 
  {
    double t = t0 + xi*dt; 
    double tn = xj[xi]; 
    x[xi] = t; 
    if (tn == t) 
    {
      y[xi] = g->GetY()[xi]; 
      continue; 
    }
    tn = (tn - t0) / dt; 

    double ans = g->GetY()[xi]; 
    ans *= (xi % 2) ? -pi : pi;
    

    int skip = -1; 
    // check if tn is another sample value... if so will cancel out part of product
    ans *= (tn -xi); 
    if (fabs(round(tn) - tn) < 1e-9 && tn >=0  && tn < n )
    {
      ans = 0; 
      int j = round(tn); 
      skip = j; 
      double tk = (xj[j] - t0)/dt; 
      ans *= (j %2) ? -1./pi : 1./pi; 
      ans *=  (xi-tk)  / (xi-j) / (tn-tk);
    }
    else
    {
      ans /= sin(pi * tn); 
    }

    for (int j = 0; j < n; j++) 
    {
      if (j == xi || j == skip) continue; 
      double tk = (xj[j] - t0)/dt; 
      ans *= (tn -j) * (xi-tk)  / (xi-j) / (tn-tk);
    }
 


    for (int i = 0; i < n; i++) 
    {
      if (i == xi) continue; 
      double dans = g->GetY()[i]; 
      dans *= (xi % 2) ? -pi : pi; 
      dans /=(xi-i); 
      double tn = (xj[i]-t0)/dt; 
      skip = -1; 

      if (fabs(tn - i) < 1e-9)
      {
          dans *= (i % 2) ? -1./pi :  1./pi; 
      }
      else if (fabs(round(tn) -tn) < 1e-9)
      {
        skip = round(tn); 
        dans *= (skip % 2) ? -1./pi :  1./pi; 
        dans *= (tn-i); 
      }
      else
      {
          dans /= sin(pi * tn) / (tn - i); 
      }

 
      for (int j = 0; j < n; j++) 
      {
        if (j == i) continue; 
        double tk = (xj[j] - t0)/dt; 
        dans *= (xi-tk)  / (tn-tk);
        if (j!=skip) dans *= (tn -j); 
        if (j!=xi) dans /= (xi -j); 
      }

      ans += dans;
    }

    y[xi] = ans; 

  }

//  w.Print(); 
  TGraph * gout = new TGraph(n,x,y); 
  if (supersample!=1) 
  {
    TGraph * grealout = FFTtools::supersample(gout,supersample); 
    delete gout; 
    return grealout;
  }
  return gout; 
}


TGraph * FFTtools::getInterpolatedGraphSincKernel(const TGraph *g, double dt, int nout)
{

  const double * xj = g->GetX(); 
  const double * yj = g->GetY(); 
  int n = g->GetN(); 
  infer_vals(dt,n,nout,xj); 
  
  double jac[n]; 
  jac[0] = (xj[n-1] - xj[0])/(n-1); 

  double yout[nout]; 
  double xout[nout]; 

  for (int i = 1; i < n; i++) 
  {
    jac[i] = xj[i] - xj[i-1]; 
  }

  for (int i = 0; i < nout; i++) 
  {
    xout[i] = xj[0] + i * dt; 
    yout[i] = 0; 

    for (int j = 0; j < n; j++) 
    {
      yout[i] += jac[j]/jac[0] * yj[j]*sinc((xout[i]-xj[j])/dt); 
    }
  }

  return new TGraph (nout,xout,yout); 
}



FFTWComplex * FFTtools::getUnevenDFT(const TGraph *g, double df, int nout)
{
  int n=g->GetN();
  const double *xj= g->GetX();
  const double *yj= g->GetY();

  FFTWComplex *dft = new FFTWComplex [nout/2+1];
  double dt =  1./ (double(nout) * df); 

  for(int k=0;k<nout/2+1;k++){
    dft[k].re=0.;
    dft[k].im=0.;
    double freq = df*k;
    for(int j=0;j<n;j++){
      double arg=2.*TMath::Pi()*freq*xj[j];
      double xlast = j == 0 ? -dt : xj[j-1]; 
      double xnext = j == n-1 ? xj[j]+dt : xj[j+1]; 
      double weight = sqrt(0.5 *(xnext-xlast)/dt); 
      dft[k].re+=yj[j]*cos(-arg)*weight;
      dft[k].im+=yj[j]*sin(-arg)*weight; 
    }
  }
  return dft;
}

TGraph *FFTtools::getInterpolatedGraphDFT(const TGraph *g, double dt, int nout, double maxF) 
{
  const double * xj = g->GetX(); 
  int n=  g->GetN();
  infer_vals(dt,n,nout,xj); 

  double df =  1./ (double(nout) * dt); 
  FFTWComplex * dft = getUnevenDFT(g, df, nout); 
  if (maxF)
  {
    int cutoff = maxF / df + 0.5; 
    for (int  i =cutoff; i < nout/2 + 1; i++) 
    {

      dft[i] = FFTWComplex(0,0); 
    }
 
  }
  double *new_y=FFTtools::doInvFFT(nout, dft);
  double *new_x = new double [nout];
  for (int i = 0; i < nout; i++) 
  {
    new_x[i] = dt * i; 
  }

  TGraph * new_g = new TGraph(nout,new_x,new_y); 
  delete new_x;;
  delete new_y; 

  return new_g;

}

