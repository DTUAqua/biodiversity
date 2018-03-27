#define TMB_LIB_INIT R_init_biodiversity
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(code)
  DATA_VECTOR(lat);
  DATA_VECTOR(lon);
  DATA_VECTOR(nsp);
  DATA_VECTOR_INDICATOR(keep, nsp);
  DATA_VECTOR(abund);
  DATA_VECTOR(asampl);
  DATA_VECTOR(temp);
  DATA_IVECTOR(cat);
  DATA_VECTOR(npp);
  DATA_VECTOR(mesh);  
  DATA_VECTOR(siz);
  DATA_VECTOR(depth);
  DATA_VECTOR(density);
  Type nll=0;
  int nobs=nsp.size();
  vector<Type> mu(nobs);
  Type k;
  if(code==(-2)){
    PARAMETER(loga7);
    PARAMETER(logk);
    Type a7=exp(loga7);
    k=exp(logk);  
    for(int i=0; i<nobs; ++i){
      mu(i)=a7;
    }
  }
  if(code==(-1)){
    PARAMETER_VECTOR(logb5);
    PARAMETER(logk);
    vector<Type> b5=exp(logb5);
    k=exp(logk);
    for(int i=0; i<nobs; ++i){
      mu(i)=b5(i)+Type(1.0e-12);
    }
  }
  if(code==1){  
    PARAMETER_VECTOR(b5);
    PARAMETER(b0);
    PARAMETER(b1);
    PARAMETER(b2);
    PARAMETER(logb3);
    PARAMETER(b4);
    PARAMETER(logk);
    vector<Type> a=exp(b0+b1*siz*siz);
    vector<Type> nu=exp(-b2/temp);
    Type b3=exp(logb3);
    k=exp(logk);
    
    for(int i=0; i<nobs; ++i){
      mu(i)=a(i)*nu(i)*
        log(Type(1)-(b3*logb3/(1-b3))*(abund(i))/(a(i)*nu(i)))*
	exp(b4*log(asampl(i)))*exp(b5(cat(i))*log(mesh(i)));
    }  
  }
  if(code==2){
    PARAMETER(loga7);
    PARAMETER_VECTOR(b8);
    PARAMETER(logb1);
    PARAMETER(logb2);
    PARAMETER(logb3);
    PARAMETER(logb5);
    PARAMETER(logb6);
    PARAMETER(logb7);
    PARAMETER(logk);
    Type b1=exp(logb1);
    Type b2=exp(logb2);
    Type b3=-exp(logb3);
    Type b5=exp(logb5);
    Type b6=exp(logb6);
    Type b7=exp(logb7);
    k=exp(logk);
    vector<Type> a=exp(loga7-b7*siz*siz);
    vector<Type> nu=exp(-b1/temp); 
    for(int i=0; i<nobs; ++i){
      mu(i)=a(i)*nu(i)*
      pow(npp(i),b2)*pow(depth(i),b3)*pow(asampl(i),b5)*pow(density(i),b6)*
	exp(b8(cat(i))*log(mesh(i)))+Type(1.0e-12);
    }
  } 
  if(code==3){    
    PARAMETER(loga7);
    PARAMETER_VECTOR(b5);
    PARAMETER(b0);
    PARAMETER(b1);
    PARAMETER(b2);
    PARAMETER(b4);
    PARAMETER(logk);

    vector<Type> loga=(siz+b1*siz*siz);
    vector<Type> a=exp(loga);
    Type a7=exp(loga7);
    k=exp(logk);
    for(int i=0; i<nobs; ++i){
      mu(i)=a7*a(i)*
        pow(lat(i), b2)*pow(lon(i), b0)*
        pow(asampl(i),b4)*exp(b5(cat(i))*log(mesh(i)))+Type(1.0e-12);
    }
  }
  if(code==4){    
    PARAMETER(loga7);
    PARAMETER_VECTOR(b5);
    PARAMETER(logb1);
    PARAMETER(logb2);
    PARAMETER(logb3);
    PARAMETER(logb4);
    PARAMETER(logk);

    Type a7=exp(loga7);
    Type b1=exp(logb1);
    Type b2=exp(logb2);
    Type b3=exp(logb3);
    Type b4=exp(logb4);
    k=exp(logk);
  
    vector<Type> nu=exp(-b1/temp);
    vector<Type> esiz=exp(siz);
 
    for(int i=0; i<nobs; ++i){
      mu(i)=a7*nu(i)*pow(density(i),b2)*pow(esiz(i),b3)*pow(asampl(i),b4)*pow(mesh(i),b5(cat(i)))+Type(1.0e-12);
    }  
  }

  vector<Type> var = mu + mu*mu/(k+Type(1.0e-12));
  vector<Type> logL= dnbinom2(nsp,mu,var,true);
  nll = -sum(vector<Type>(keep*logL));

  REPORT(mu);
  REPORT(var);
  REPORT(nll);
  REPORT(logL);
  return nll;
}
