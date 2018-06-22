#include <Rcpp.h>
using namespace Rcpp;


//' Extracting factors for conditional log-linear model
//' 
//' This is an auxiliary function to print out the factors for conditional log-linear model
//' given edge information.
//' 
//' 
//' @param pars a set of parameters
//' @param newcombined a \code{[(2+nc) x m ]} matrix comprised of outcomes (first row), treatments (second row), and confounders (from the third row), where \code{nc} is the number of confounders.
//' @param edgeY a matrix of which each row indicates a pair of index for adjacent outcomes.
//' @param edgeAY a matrix of which each row indicates a index for treatment (first column) and for outcome (second column) on which the treatment has a direct effect.
//' @param edgeExtra a list of edges of which a list of matrix specifying additional directed edges (from confounders or treatment to the outcomes) information.
//' 
//' 
//' @return a sum of factors.
//' 
//' @useDynLib netchain
//' @importFrom Rcpp sourceCpp
//' 
//' @export
// [[Rcpp::export]]
double multimainfunction(NumericVector pars, NumericMatrix newcombined,
                         NumericMatrix edgeY,
                         NumericMatrix edgeAY, List edgeExtra){
  

  IntegerVector lengthY = seq(1, newcombined.ncol());
  NumericVector Y = newcombined.row(0);
  NumericVector A = newcombined.row(1);
  double mainpart = sum(Y[lengthY-1]*pars[lengthY-1]);
  for(int i = 0; i < edgeY.nrow() ; ++i){
    mainpart += Y[edgeY(i,0)-1]*Y[edgeY(i,1)-1]*pars[newcombined.ncol()+i];
  }
  for(int j=0; j < edgeAY.nrow(); ++j){
    mainpart += A[edgeAY(j,0)-1]*Y[edgeAY(j,1)-1]*pars[newcombined.ncol()+edgeY.nrow()+j];
  }
  
  for(int k=0; k < edgeExtra.size(); ++k){
    IntegerMatrix tmp = edgeExtra[k];
    double factors = 1;
    for(int ll=0; ll < tmp.nrow(); ll++){
      factors *= newcombined(tmp(ll,0)-1, tmp(ll,1)-1);
    }
    mainpart += factors*pars[newcombined.ncol()+edgeY.nrow()+edgeAY.nrow()+k];
  }
  
  return(mainpart);
}


//' Calculating normalizing constant in conditional log-linear model.
//'
//' @param pars a set of parameters
//' @param combined a \code{[(2+nc) x m ]} matrix comprised of outcomes (first row), treatments (second row), and confounders (from the third row), where \code{nc} is the number of confounders.
//' @param permutetab a matrix comprised of every possible values for outcome in each row. 
//' @param edgeY a matrix of which each row indicates a pair of index for adjacent outcomes.
//' @param edgeAY a matrix of which each row indicates a index for treatment (first column) and for outcome (second column) on which the treatment has a direct effect.
//' @param edgeExtra a list of edges of which a list of matrix specifying additional directed edges (from confounders or treatment to the outcomes) information.
//' 
//' @return a normalizing constant
//' @export
//' @useDynLib netchain
//' @importFrom Rcpp sourceCpp
//'
//[[Rcpp::export]]
double multipartition(NumericVector pars, NumericMatrix combined,
                      NumericMatrix permutetab, NumericMatrix edgeY,
                      NumericMatrix edgeAY, List edgeExtra){
  
  
  double part = 0;
  for(int p=0; p < permutetab.nrow(); ++p){
    NumericMatrix dummyobs = clone(combined);
    dummyobs.row(0) = permutetab.row(p);
    part += exp(multimainfunction(pars, dummyobs, edgeY, edgeAY, edgeExtra));
  }
  return(part);
}


//' Derive log-likelihood of conditional log-linear model given parameters.
//'
//' @param pars a set of parameters
//' @param listobservations a collection of \code{[(2+nc) x m ]} matrices comprised of outcomes (first row), treatments (second row), and confounders (from the third row), where \code{nc} is the number of confounders. 
//' @param permutetab a matrix comprised of every possible values for outcome in each row. 
//' @param edgeY a matrix of which each row indicates a pair of index for adjacent outcomes.
//' @param edgeAY a matrix of which each row indicates a index for treatment (first column) and for outcome (second column) on which the treatment has a direct effect.
//' @param edgeExtra a list of edges of which a list of matrix specifying additional directed edges (from confounders or treatment to the outcomes) information.
//'
//' @return log-likelihood of conditional log-linear model given parameters, observations, and edge information.
//' @export
//'
//' @useDynLib netchain
//' @importFrom Rcpp sourceCpp
//'
//[[Rcpp::export]]
double multiloglikechain(NumericVector pars, List listobservations,
                         NumericMatrix permutetab, NumericMatrix edgeY,
                         NumericMatrix edgeAY, List edgeExtra){
  
  
  double whole = 0;
  double Zpart = 0;
  
  for(int ii=0; ii < listobservations.size(); ++ii){
    NumericMatrix tmpmatrix = listobservations[ii];
    whole += multimainfunction(pars, tmpmatrix, edgeY, edgeAY, edgeExtra);
    Zpart += log(multipartition(pars, tmpmatrix, permutetab, edgeY, edgeAY, edgeExtra));
  }
  
  double loglike = whole - Zpart;
  return(loglike);
}


