#include "DPWignerFunctionGeneral.hh"
#include "TMath.h"
#include <iostream>

 long DPWignerFunctionGeneral::factorials[]={1,
1,
2,
6,
24,
120,
720,
5040,
40320,
362880,
3628800,
39916800,
479001600,
6227020800,
87178291200,
1307674368000,
20922789888000,
355687428096000,
6402373705728000,
121645100408832000};

long DPWignerFunctionGeneral::fact(int n)
{
  // Up to factorial 19 we have precalculated values
  if ( n<=19 )
  {
    return factorials[n];
  }
  else
  {
    return n*fact(n-1);
  }
}

double DPWignerFunctionGeneral::function(double cosTheta, double n, double m)
{
  double result=0;

  int jmin=0;
  if (m-n>jmin)
  {
    jmin=m-n;
  }
  int jmax=j-n;
  if (j+m<jmax)
  {
    jmax=j+m;
  }

  double prefactor=fact(j+n)*fact(j-n)*fact(j+m)*fact(j-m);
  prefactor=TMath::Sqrt(prefactor);

//  std::cout<<"Wigner: "<<j<<" "<<n<<" "<<m<<" ";
//  std::cout<<jmin<<" "<<jmax<<" "<<prefactor<<std::endl;

  double theta=TMath::ACos(cosTheta);
  double cos=TMath::Cos(theta/2.);
  double sin=TMath::Sin(theta/2.);

  for (int s=jmin;s<=jmax;++s)
  {
    result+=TMath::Power(-1,n-m+s)/
       (fact(j+m-s)*fact(s)*fact(n-m+s)*fact(j-n-s))*
       TMath::Power(cos,2*j+m-n-2*s)*
       TMath::Power(sin,n-m+2*s);
//    std::cout<<"Wigner s-loop: "<<s<<" "<<n-m+s<<" "<<(fact(j+m-s)*fact(s)*fact(n-m+s)*fact(j-n-s))
//             <<" "<<2*j+m-n-2*s<<" "<<n-m+2*s<<std::endl;
  }
  return prefactor*result;
}

