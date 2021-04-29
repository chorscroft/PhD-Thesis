


statsFormula<-c(
  expression(paste("Z",''[alpha])),
  expression(paste("Z",''[beta])),
  expression(paste(bgroup("(",atop("|L|",2),")"),"+",bgroup("(",atop("|R|",2),")"))),
  "|L||R|",
  expression(paste("Z",''[alpha],"+","Z",''[beta])), #"Zalpha_plus_Zbeta",
  expression(paste("Z",''[alpha],"-","Z",''[beta])), #"Zalpha_minus_Zbeta",
  expression(paste("Z",''[alpha],"Z",''[beta])), #"ZalphaZbeta",
  expression(paste(frac(paste("Z",''[alpha]),paste("Z",''[beta])))), #"Zalpha_over_Zbeta",
  expression(paste("Z"[alpha]^paste(r^2,"/E[",r^2,"]"))), #"Zalpha_rsq_over_expected",
  expression(paste("Z"[alpha]^paste("log(",r^2,"/E[",r^2,"])"))), #"Zalpha_log_rsq_over_expected",
  expression(paste("Z"[alpha]^ZScore)), #"Zalpha_Zscore",
  expression(paste("Z"[alpha]^BetaCDF)), #"Zalpha_BetaCDF",
  expression(paste("Z"[alpha],"-","Z"[alpha]^paste("E[",r^2,"]"))), #"Zalpha_minus_Zalpha_expected",
  expression(paste(frac("Z"[alpha],"Z"[alpha]^paste("E[",r^2,"]")))), #"Zalpha_over_Zalpha_expected",
  expression(paste("Z"[beta]^paste(r^2,"/E[",r^2,"]"))), #"Zbeta_rsq_over_expected",
  expression(paste("Z"[beta]^paste("log(",r^2,"/E[",r^2,"])"))), #"Zbeta_log_rsq_over_expected",
  expression(paste("Z"[beta]^ZScore)), #"Zbeta_Zscore",
  expression(paste("Z"[beta]^BetaCDF)), #"Zbeta_BetaCDF",
  expression(paste("Z"[beta],"-","Z"[beta]^paste("E[",r^2,"]"))), #"Zbeta_minus_Zbeta_expected",
  expression(paste(frac("Z"[beta],"Z"[beta]^paste("E[",r^2,"]")))), #"Zbeta_over_Zbeta_expected",
  expression(paste("Z"[alpha]^paste(r^2,"/E[",r^2,"]"),"+","Z"[beta]^paste(r^2,"/E[",r^2,"]"))), #"Zalpha_rsq_over_expected_plus_Zbeta_rsq_over_expected",
  expression(paste("Z"[alpha]^paste("log(",r^2,"/E[",r^2,"])"),"+","Z"[beta]^paste("log(",r^2,"/E[",r^2,"])"))), #"Zalpha_log_rsq_over_expected_plus_Zbeta_log_rsq_over_expected",
  expression(paste("Z"[alpha]^ZScore,"+","Z"[beta]^ZScore)), #"Zalpha_Zscore_plus_Zbeta_Zscore",
  expression(paste("Z"[alpha]^BetaCDF,"+","Z"[beta]^BetaCDF)), #"Zalpha_BetaCDF_plus_Zbeta_BetaCDF",
  expression(paste("(Z"[alpha],"-","Z"[alpha]^paste("E[",r^2,"]"),")+","(Z"[beta],"-","Z"[beta]^paste("E[",r^2,"]"),")")), #"Zalpha_minus_Zalpha_expected_plus_Zbeta_minus_Zbeta_expected",
  expression(paste(frac("Z"[alpha],"Z"[alpha]^paste("E[",r^2,"]")),"+",frac("Z"[beta],"Z"[beta]^paste("E[",r^2,"]")))), #"Zalpha_over_Zalpha_expected_plus_Zbeta_over_Zbeta_expected",
  expression(paste("Z"[alpha]^paste(r^2,"/E[",r^2,"]"),"-","Z"[beta]^paste(r^2,"/E[",r^2,"]"))), #"Zalpha_rsq_over_expected_minus_Zbeta_rsq_over_expected",
  expression(paste("Z"[alpha]^paste("log(",r^2,"/E[",r^2,"])"),"-","Z"[beta]^paste("log(",r^2,"/E[",r^2,"])"))), #"Zalpha_log_rsq_over_expected_minus_Zbeta_log_rsq_over_expected",
  expression(paste("Z"[alpha]^ZScore,"-","Z"[beta]^ZScore)), #"Zalpha_Zscore_minus_Zbeta_Zscore",
  expression(paste("Z"[alpha]^BetaCDF,"-","Z"[beta]^BetaCDF)), #"Zalpha_BetaCDF_minus_Zbeta_BetaCDF",
  expression(paste("(Z"[alpha],"-","Z"[alpha]^paste("E[",r^2,"]"),")-","(Z"[beta],"-","Z"[beta]^paste("E[",r^2,"]"),")")), #"Zalpha_minus_Zalpha_expected_minus_Zbeta_minus_Zbeta_expected",
  expression(paste(frac("Z"[alpha],"Z"[alpha]^paste("E[",r^2,"]")),"-",frac("Z"[beta],"Z"[beta]^paste("E[",r^2,"]")))), #"Zalpha_over_Zalpha_expected_minus_Zbeta_over_Zbeta_expected",
  expression(paste("Z"[alpha]^BetaCDF,"Z"[beta]^BetaCDF)), #"Zalpha_BetaCDFZbeta_BetaCDF",
  expression(paste(frac("Z"[alpha]^BetaCDF,"Z"[beta]^BetaCDF))), #"Zalpha_BetaCDF_over_Zbeta_BetaCDF",
  expression(paste(frac(paste("Z",''[alpha]),paste("Z",''[beta])),"-",frac("Z"[alpha]^paste("E[",r^2,"]"),"Z"[beta]^paste("E[",r^2,"]")))) #"Zalpha_over_Zbeta_minus_Zalpha_expected_over_Zbeta_expected"
)

statsList<-list(
  Zalpha=expression(paste("Z",''[alpha])),
  Zbeta=expression(paste("Z",''[beta])),
  L_plus_R=expression(paste(bgroup("(",atop("|L|",2),")"),"+",bgroup("(",atop("|R|",2),")"))),
  LR="|L||R|",
  Zalpha_plus_Zbeta=expression(paste("Z",''[alpha],"+","Z",''[beta])),
  Zalpha_minus_Zbeta=expression(paste("Z",''[alpha],"-","Z",''[beta])),
  ZalphaZbeta=expression(paste("Z",''[alpha],"Z",''[beta])),
  Zalpha_over_Zbeta=expression(paste(frac(paste("Z",''[alpha]),paste("Z",''[beta])))),
  Zalpha_rsq_over_expected=expression(paste("Z"[alpha]^paste(r^2,"/E[",r^2,"]"))),
  Zalpha_log_rsq_over_expected=expression(paste("Z"[alpha]^paste("log(",r^2,"/E[",r^2,"])"))),
  Zalpha_Zscore=expression(paste("Z"[alpha]^ZScore)),
  Zalpha_BetaCDF=expression(paste("Z"[alpha]^BetaCDF)),
  Zalpha_minus_Zalpha_expected=expression(paste("Z"[alpha],"-","Z"[alpha]^paste("E[",r^2,"]"))),
  Zalpha_over_Zalpha_expected=expression(paste(frac("Z"[alpha],"Z"[alpha]^paste("E[",r^2,"]")))),
  Zbeta_rsq_over_expected=expression(paste("Z"[beta]^paste(r^2,"/E[",r^2,"]"))),
  Zbeta_log_rsq_over_expected=expression(paste("Z"[beta]^paste("log(",r^2,"/E[",r^2,"])"))),
  Zbeta_Zscore=expression(paste("Z"[beta]^ZScore)),
  Zbeta_BetaCDF=expression(paste("Z"[beta]^BetaCDF)),
  Zbeta_minus_Zbeta_expected=expression(paste("Z"[beta],"-","Z"[beta]^paste("E[",r^2,"]"))),
  Zbeta_over_Zbeta_expected=expression(paste(frac("Z"[beta],"Z"[beta]^paste("E[",r^2,"]")))),
  Zalpha_rsq_over_expected_plus_Zbeta_rsq_over_expected=expression(paste("Z"[alpha]^paste(r^2,"/E[",r^2,"]"),"+","Z"[beta]^paste(r^2,"/E[",r^2,"]"))),
  Zalpha_log_rsq_over_expected_plus_Zbeta_log_rsq_over_expected=expression(paste("Z"[alpha]^paste("log(",r^2,"/E[",r^2,"])"),"+","Z"[beta]^paste("log(",r^2,"/E[",r^2,"])"))),
  Zalpha_Zscore_plus_Zbeta_Zscore=expression(paste("Z"[alpha]^ZScore,"+","Z"[beta]^ZScore)),
  Zalpha_BetaCDF_plus_Zbeta_BetaCDF=expression(paste("Z"[alpha]^BetaCDF,"+","Z"[beta]^BetaCDF)),
  Zalpha_minus_Zalpha_expected_plus_Zbeta_minus_Zbeta_expected=expression(paste("(Z"[alpha],"-","Z"[alpha]^paste("E[",r^2,"]"),")+","(Z"[beta],"-","Z"[beta]^paste("E[",r^2,"]"),")")),
  Zalpha_over_Zalpha_expected_plus_Zbeta_over_Zbeta_expected=expression(paste(frac("Z"[alpha],"Z"[alpha]^paste("E[",r^2,"]")),"+",frac("Z"[beta],"Z"[beta]^paste("E[",r^2,"]")))),
  Zalpha_rsq_over_expected_minus_Zbeta_rsq_over_expected=expression(paste("Z"[alpha]^paste(r^2,"/E[",r^2,"]"),"-","Z"[beta]^paste(r^2,"/E[",r^2,"]"))),
  Zalpha_log_rsq_over_expected_minus_Zbeta_log_rsq_over_expected=expression(paste("Z"[alpha]^paste("log(",r^2,"/E[",r^2,"])"),"-","Z"[beta]^paste("log(",r^2,"/E[",r^2,"])"))),
  Zalpha_Zscore_minus_Zbeta_Zscore=expression(paste("Z"[alpha]^ZScore,"-","Z"[beta]^ZScore)),
  Zalpha_BetaCDF_minus_Zbeta_BetaCDF=expression(paste("Z"[alpha]^BetaCDF,"-","Z"[beta]^BetaCDF)),
  Zalpha_minus_Zalpha_expected_minus_Zbeta_minus_Zbeta_expected=expression(paste("(Z"[alpha],"-","Z"[alpha]^paste("E[",r^2,"]"),")-","(Z"[beta],"-","Z"[beta]^paste("E[",r^2,"]"),")")),
  Zalpha_over_Zalpha_expected_minus_Zbeta_over_Zbeta_expected=expression(paste(frac("Z"[alpha],"Z"[alpha]^paste("E[",r^2,"]")),"-",frac("Z"[beta],"Z"[beta]^paste("E[",r^2,"]")))),
  Zalpha_BetaCDFZbeta_BetaCDF=expression(paste("Z"[alpha]^BetaCDF,"Z"[beta]^BetaCDF)),
  Zalpha_BetaCDF_over_Zbeta_BetaCDF=expression(paste(frac("Z"[alpha]^BetaCDF,"Z"[beta]^BetaCDF))),
  Zalpha_over_Zbeta_minus_Zalpha_expected_over_Zbeta_expected=expression(paste(frac(paste("Z",''[alpha]),paste("Z",''[beta])),"-",frac("Z"[alpha]^paste("E[",r^2,"]"),"Z"[beta]^paste("E[",r^2,"]"))))
)





plot(1,2,xlab=statsFormula[35])
