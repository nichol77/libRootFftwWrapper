

void testAkima() {
  //http://www.iue.tuwien.ac.at/phd/rottinger/node60.html
  Double_t a[4];

  Double_t dx=(thisX-x[i]);
  Double_t s =a[0] + a[1] * dx + a[2]*dx*dx + a[3]*dx*dx*dx; 
  Double_t sdashx = a[1] + 2*a[2]*dx + 3*a[3]*dx*dx; 
  
  
  Double_t si=a[0]; //Since dxi =0, so s(x_i)=a[0]  this is y[i]
  Double_t Dx=(x[i+1]-x[i]);
  Double_t sip1=a[0] + a[1] * Dx + a[2]*Dx*Dx + a[3]*Dx*Dx*Dx; 
  Double_t sdashi=a[1];  
  Double_t sdaship1=a[1] + 2*a[2]*Dx + 3*a[3]*Dx*Dx; 


  Double_t dydx[N+3] = (y[i+1] - y[i]) / (x[i+1] - x[i]);
  dydx[0]=(2*dydx[1])-dydx[2];
  dydx[-1]=(2*dydx[0])-dydx[1]; //Will stick this in N+2
  dydx[N]=(2*dydx[N-1])-dydx[N-2];
  dydx[N+1]=(2*dydx[N])-dydx[N-1];

  Double_t raww[N+3]= TMath::Abs(dydx[i] - dydx[i-1]);
  
  Double_t w_i = raww[i-1];
  Double_t w_im1 = raww[i+1];

  Double_t sdashi = (w_im1 * dydx[i-1] + w_i* dydx[i]) / ( w_im1 +w_i)
  
  
}
