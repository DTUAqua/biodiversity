#define TMB_LIB_INIT R_init_biodiversity
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_VECTOR(lat);
  DATA_VECTOR(nsp);
  DATA_VECTOR_INDICATOR(keep, nsp);
  DATA_VECTOR(abund);
  DATA_VECTOR(logsiz);
  DATA_VECTOR(asampl);
  DATA_VECTOR(temp);
  DATA_IVECTOR(cat);
  DATA_VECTOR(npp);
  DATA_VECTOR(mesh);  
  
  PARAMETER_VECTOR(b5);
  PARAMETER(b0);
  PARAMETER(b1);
  PARAMETER(b2);
  PARAMETER(logb3);
  PARAMETER(b4);
  PARAMETER(logk);
  vector<Type> a=exp(b0+b1*logsiz*logsiz);
  vector<Type> nu=exp(-b2/temp);
  Type b3=exp(logb3);
  Type k=exp(logk);
  
  int nobs=nsp.size();
  vector<Type> mu(nobs);
  
  for(int i=0; i<nobs; ++i){
      mu(i)=a(i)*nu(i)*
        log(Type(1)-(b3*logb3/(1-b3))*(abund(i))/(a(i)*nu(i)))*
	exp(b4*log(asampl(i)))*exp(b5(cat(i))*log(mesh(i)));
      }
  
  vector<Type> var = mu + mu*mu/k;
  Type nll=0;
  for(int i=0; i<nobs; ++i){
    nll += -keep(i)*dnbinom2(nsp(i),mu(i),var(i),true);
  }
  
  REPORT(mu);
  REPORT(var);
  REPORT(b0);
  REPORT(b1);
  REPORT(b2);
  REPORT(b3);
  REPORT(b4);
  REPORT(b5);
  REPORT(k);
    
  ADREPORT(b0);
  ADREPORT(b1);
  ADREPORT(b2);
  ADREPORT(b3);
  ADREPORT(b4);
  ADREPORT(b5);
  ADREPORT(k);
  REPORT(nll);
  return nll;
}

