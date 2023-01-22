# p0m: probability of misleading evidence under the null model 
# p1m: probability of misleading evidence under the alternative model 
# pc: overall probability of correct evidence 
# pi0: prior probability of the null model
# mm: number of sites m 
# S: size of the random sample generated from the analysis prior 
# T: size of the random sample from the predictive distribution of BF_{01} 
# s_a: scale parameter of the analysis prior
# mu_d: location parameter of the design prior
# s_d: scale parameter of the design prior

### a function that given m, p0m, p1m, pc, pi0 and the prior parameters ###
### returns the minimum n that guarantees the predetermined pc ### 


n_uncond=function(p0m=0.05, p1m=0.05, pc=0.8, pi0=0.5, mm=8, S=10^4, T=5*10^4,  s_a=1/7, mu_d=0.2, s_d=1/55 ){
  
  RESULTS=data.frame(n=NA, m=mm, p1c=NA, rk1=NA, k0=NA)
  
  BF=function(Qm, n, m, t) {
    
    bfa=c()
    for (i in 1:T) {
      
      BFa=(n^((m-1)/2))*exp(-Qm[i]/2)/mean( (1/n+t^2)^((1-m)/2) * exp(-Qm[i]/2*(1+n*t^2)^(-1) ) )
      bfa=c(bfa,BFa)
      
    }  
    
    return(bfa)
    
  }
  
  set.seed(1)
  a0=rt(S, df=4)
  aa2=abs(s_a*a0) 
  
  set.seed(2)
  a=rt(T, df=4)
  aa=abs(s_d*a+mu_d) 
  
  ind0=0
  for (mi in mm) {
    
    ind0=ind0+1
    stop=0
    ## first initial value (n)
    ind=1
    ni1=round(500*pc^2/mi)
    

    set.seed(3)
    nulldata=rchisq(T, df=mi-1)
    
    set.seed(4)
    idata=rchisq(T, df=mi-1)
    altdata=idata*(1+ni1*aa^2)
    
    
    bf0=BF(nulldata, n=ni1, m=mi, t=aa2)
    bf1=BF(altdata, n=ni1, m=mi, t=aa2)
    
    k0=quantile(bf1, 1-p1m)
    k1=quantile(bf0, p0m)
    
    p1c=sum(bf1<(k1)) / length(bf1)
    p0c = sum(bf0 > k0) / length(bf0) 
    
    pc1=pi0*p0c+(1-pi0)*p1c
    
    print(paste0("It ",ind,": Suggested (n,m)=(", ni1,",",mi,") and achieved pc=",pc1))
    
    if (abs(pc1-pc)/pc<0.001) {stop=1 ; ninew=ni1 ; pcnew=pc1}
    
    
    if (stop == 0)  {  
      ## second initial value (n)
      ind=2
      ni2=2*ni1
      
      
      altdata=idata*(1+ni2*aa^2)
      
      
      bf0=BF(nulldata, n=ni2, m=mi, t=aa2)
      bf1=BF(altdata, n=ni2, m=mi, t=aa2)
      
      k0=quantile(bf1, 1-p1m)
      k1=quantile(bf0, p0m)
      
      p1c =sum(bf1<(k1)) / length(bf1)
      p0c = sum(bf0 > k0) / length(bf0) 
      
      pc2=pi0*p0c+(1-pi0)*p1c
      
      print(paste0("It ",ind,": Suggested (n,m)=(", ni2,",",mi,") and achieved pc=",pc2))
      
      if (abs(pc2-pc)/pc<0.001) {stop=1 ; ninew=ni2 ; pcnew=pc2}
      
      
      #########################################################################################
      
      while(pc2==pc1) {
        ni2=ni2+1
        ind=ind+1
        
        
        altdata=idata*(1+ni2*aa^2)
        
        
        bf0=BF(nulldata, n=ni2, m=mi, t=aa2)
        bf1=BF(altdata, n=ni2, m=mi, t=aa2)
        
        k0=quantile(bf1, 1-p1m)
        k1=quantile(bf0, p0m)
        
        p1c=sum(bf1<(k1)) / length(bf1)
        p0c = sum(bf0 > k0) / length(bf0) 
        
        pc2=(p0c+pc1)/2
        
        print(paste0("It ",ind,": Suggested (n,m)=(", ni2,",",mi,") and achieved pc=",pc2))
      }
      
      if (abs(pc2-pc)/pc<0.001) {stop=1 ; ninew=ni2 ; pcnew=pc2}
      
      
    }
    #########################################################################################
    
    if (stop==0) {
      
      ninew=max(round(ni2+(pc-pc2)*(ni2-ni1)/(pc2-pc1)),1)    
      ind=ind+1
      
      
      altdata=idata*(1+ninew*aa^2)
      
      
      bf0=BF(nulldata, n=ninew, m=mi, t=aa2)
      bf1=BF(altdata, n=ninew, m=mi, t=aa2)
      
      k0=quantile(bf1, 1-p1m)
      k1=quantile(bf0, p0m)
      
      p1c=sum(bf1<(k1)) / length(bf1)
      p0c = sum(bf0 > k0) / length(bf0) 
      
      pcnew=pi0*p0c+(1-pi0)*p1c
      
      print(paste0("It ",ind,": Suggested (n,m)=(", ninew,",",mi,") and achieved pc=",pcnew))
      
      vn=c(ni1,ni2,ninew) 
      if (length(unique(vn)) == 1) { 
        message("Warning message: the desired accuracy cannot be achieved.")
        return(hnew)
      }
      
      while( abs(pcnew-pc)/pc>0.001) {
        ni1=ni2
        pc1=pc2
        ni2=ninew
        pc2=pcnew
        
        while(pc2==pc1) {
          ni2=ni2+1
          ind=ind+1
          
          altdata=idata*(1+ni2*aa^2)
          
          
          bf0=BF(nulldata, n=ni2, m=mi, t=aa2)
          bf1=BF(altdata, n=ni2, m=mi, t=aa2)
          
          k0=quantile(bf1, 1-p1m)
          k1=quantile(bf0, p0m)
          
          p1c=sum(bf1<(k1)) / length(bf1)
          p0c = sum(bf0 > k0) / length(bf0) 
          
          pc2=pi0*p0c+(1-pi0)*p1c
          
          print(paste0("It ",ind,": Suggested (n,m)=(", ni2,",",mi,") and achieved pc=",pc2))
        }
        
        
        ninew=max(round(ni2+(pc-pc2)*(ni2-ni1)/(pc2-pc1)),1)    
        ind=ind+1
        
        altdata=idata*(1+ninew*aa^2)
        
        
        bf0=BF(nulldata, n=ninew, m=mi, t=aa2)
        bf1=BF(altdata, n=ninew, m=mi, t=aa2)
        
        k0=quantile(bf1, 1-p1m)
        k1=quantile(bf0, p0m)
        
        p1c=sum(bf1<(k1)) / length(bf1)
        p0c = sum(bf0 > k0) / length(bf0) 
        
        pcnew=pi0*p0c+(1-pi0)*p1c
        
        print(paste0("It ",ind,": Suggested (n,m)=(", ninew,",",mi,") and achieved pc=",pcnew))
        
        
        vn=c(ni1,ni2,ninew) 
        if (length(unique(vn)) < 3) { 
          message("Warning message: the desired accuracy in Regula Falsi cannot be achieved.")
          break
        }
        
        
        
      }
    }
    
    print("End of Regula Falsi")
    #########################################################################################
    
    if (pcnew < pc) {
      
      while (pcnew < pc) {
        
        ninew=ninew+1
        ind=ind+1
        
        altdata=idata*(1+ninew*aa^2)
        
        
        bf0=BF(nulldata, n=ninew, m=mi, t=aa2)
        bf1=BF(altdata, n=ninew, m=mi, t=aa2)
        
        k0=quantile(bf1, 1-p1m)
        k1=quantile(bf0, p0m)
        
        p1c=sum(bf1<(k1)) / length(bf1)
        p0c = sum(bf0 > k0) / length(bf0) 
        
        pcnew=pi0*p0c+(1-pi0)*p1c
        
        print(paste0("It ",ind,": Suggested (n,m)=(", ninew,",",mi,") and achieved pc=",pcnew))
        
        
      } } else {
        
        while (pcnew >= pc) {
          
          niprev=ninew
          pcprev=pcnew
          
          ninew=ninew-1
          ind=ind+1
          
          altdata=idata*(1+ninew*aa^2)
          
          
          bf0=BF(nulldata, n=ninew, m=mi, t=aa2)
          bf1=BF(altdata, n=ninew, m=mi, t=aa2)
          
          k0=quantile(bf1, 1-p1m)
          k1=quantile(bf0, p0m)
          
          p1c=sum(bf1<(k1)) / length(bf1)
          p0c = sum(bf0 > k0) / length(bf0) 
          
          pcnew=pi0*p0c+(1-pi0)*p1c
          
          print(paste0("It ",ind,": Suggested (n,m)=(", ninew,",",mi,") and achieved pc=",pcnew))
          
        } 
        
        ninew=niprev
        pcnew=pcprev
        
      }
    
    RESULTS[ind0,1]=ninew
    RESULTS[ind0,3]=pcnew
    RESULTS[ind0,4]=round(k1,3)
    RESULTS[ind0,5]=round(k0,3)
    
    print(RESULTS[1:ind0,])
    
  }
  
  return(RESULTS)
  
}


n_uncond()









