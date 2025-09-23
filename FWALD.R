##############################################################################
# WALD TEST
##############################################################################

FWALD <- function(model,var){
  ind <- grep(var,names(coef(model)));
  coef <- coef(model)[ind];
  vcov <- vcov(model)[ind,ind];
  waldstat <- coef%*%solve(vcov)%*%coef;
  df <- length(coef);
  return(1-pchisq(waldstat,df));
}
