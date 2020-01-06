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
  DATA_VECTOR(totCatch);
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
      mu(i)=b5(i);//+Type(1.0e-12);
    }
  }
  if(code==1){  
    PARAMETER_VECTOR(b5);
    PARAMETER(b0);
    PARAMETER(b1);
    PARAMETER(b2);
    PARAMETER(b3);
    PARAMETER(b4);
    PARAMETER_VECTOR(b8);    
    PARAMETER(logk);
    vector<Type> a=exp(b0+b1*siz*siz);
    vector<Type> nu=exp(-b2/temp);
    Type logb3=log(b3);
    k=exp(logk);
    for(int i=0; i<nobs; ++i){
      mu(i)=a(i)*nu(i)*
        log(Type(1)-(b3*logb3/(1-b3))*(abund(i)+1.0e-12)/(a(i)*nu(i)))*
	pow(totCatch(i),b8(0))*exp(b4*log(asampl(i)))*exp(b5(cat(i))*log(mesh(i)))+Type(1.0e-12);
    }  
  }
  if(code==2){
    PARAMETER(loga7);
    PARAMETER_VECTOR(b8);
    PARAMETER(b1);
    PARAMETER(b2);
    PARAMETER(b3);
    PARAMETER(b4);    
    PARAMETER(b5);
    PARAMETER(b6);
    PARAMETER(b7);
    PARAMETER(logk);
    k=exp(logk);
    vector<Type> a=exp(loga7+b7*siz*siz);
    vector<Type> nu=exp(-b1/temp);

    for(int i=0; i<nobs; ++i){
      mu(i)=a(i)*nu(i)*
        pow(npp(i),b2)*pow(depth(i),b3)*pow(asampl(i),b5)*pow(abund(i),b6)*pow(totCatch(i),b4)*
	exp(b8(cat(i))*log(mesh(i)))+Type(1.0e-12);
    }
  } 
  if(code==3){    
    PARAMETER(loga7);
    PARAMETER_VECTOR(b5);
    PARAMETER(b0);
    PARAMETER(b1);
    PARAMETER(b2);
    PARAMETER(b3);
    PARAMETER(b4);
    PARAMETER_VECTOR(b8);    
    PARAMETER(logk);

    vector<Type> loga=(siz+b1*siz*siz);
    vector<Type> a=exp(loga);
    Type a7=exp(loga7);
    k=exp(logk);
    
    for(int i=0; i<nobs; ++i){
      mu(i)=a7*a(i)*
        pow(lat(i), b2)*pow(lon(i), b0)*pow(depth(i),b3)*
        pow(totCatch(i),b8(0))*pow(asampl(i),b4)*exp(b5(cat(i))*log(mesh(i)))+Type(1.0e-12);
    }
  }
  if(code==4){    
    PARAMETER(loga7);
    PARAMETER_VECTOR(b5);
    PARAMETER(b1);
    PARAMETER(b2);
    PARAMETER(b3);
    PARAMETER(b4);
    PARAMETER_VECTOR(b8);    
    PARAMETER(logk);

    Type a7=exp(loga7);
    k=exp(logk);
  
    vector<Type> nu=exp(-b1/temp);
    vector<Type> esiz=exp(siz);

    for(int i=0; i<nobs; ++i){
      mu(i)=a7*nu(i)*pow(density(i),b2)*pow(esiz(i),b3)*pow(totCatch(i),b8(0))*pow(asampl(i),b4)*pow(mesh(i),b5(cat(i)))+Type(1.0e-12);
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
