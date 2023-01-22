# p0m: probability of misleading evidence under the null model 
# p1c: probability of correct evidence under the alternative model 
# mm: number of sites m 
# S: size of the random sample generated from the analysis prior 
# T: size of the random sample from the predictive distribution of BF_{01} 
# s_a: scale parameter of the analysis prior
# mu_d: location parameter of the design prior
# s_d: scale parameter of the design prior

### a function that given m, p0m, p1c, and the prior parameters ###
### returns the minimum n that guarantees the predetermined p1c ### 


n_cond=function(p0m=0.05, p1c=0.8, mm=8, S=10^4, T=5*10^4,  s_a=1/7, mu_d=0.2, s_d=1/55 ){
  
  RESULTS=data.frame(n=NA, m=mm, p1c=NA, rk1=NA)
  
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
    ind=1
    ni1=round(500*p1c^2/mi)
    
    set.seed(3)
    nulldata=rchisq(T, df=mi-1)
    
    set.seed(4)
    idata=rchisq(T, df=mi-1)
    altdata=idata*(1+ni1*aa^2)
    
    
    bf0=BF(nulldata, n=ni1, m=mi, t=aa2)
    bf1=BF(altdata, n=ni1, m=mi, t=aa2)
    
    k1=quantile(bf0, p0m)
    
    p1c1=sum(bf1<(k1)) / length(bf1)

    
    print(paste0("It ",ind,": Suggested (n,m)=(", ni1,",",mi,") and achieved p1c=",p1c1))
    
    if (abs(p1c1-p1c)/p1c<0.001) {stop=1 ; ninew=ni1 ; p1cnew=p1c1}
    
    
    if (stop == 0)  {  

      ind=2
      ni2=2*ni1
      
      altdata=idata*(1+ni2*aa^2)
      
      
    bf0=BF(nulldata, n=ni1, m=mi, t=aa2)  
    bf1=BF(altdata, n=ni1, m=mi, t=aa2)
    
    k1=quantile(bf0, p0m)
    
    p1c2=sum(bf1<(k1)) / length(bf1)

  
      print(paste0("It ",ind,": Suggested (n,m)=(", ni2,",",mi,") and achieved p1c=",p1c2))
            
      if (abs(p1c2-p1c)/p1c<0.001) {stop=1 ; ninew=ni2 ; pd1cnew=p1c2}
      
      
      #########################################################################################
      
      while(p1c2==p1c1) {
        ni2=ni2+1
        ind=ind+1
        
        altdata=idata*(1+ni2*aa^2)
        
        
    bf0=BF(nulldata, n=ni1, m=mi, t=aa2)  
    bf1=BF(altdata, n=ni1, m=mi, t=aa2)
    
    k1=quantile(bf0, p0m)
    
    p1c2=sum(bf1<(k1)) / length(bf1)
        
        print(paste0("It ",ind,": Suggested (n,m)=(", ni2,",",mi,") and achieved p1c=",p1c2))
      }
      
      if (abs(p1c2-p1c)/p1c<0.001) {stop=1 ; ninew=ni2 ; p1cnew=p1c2}
      
      
    }
    #########################################################################################
    
    if (stop==0) {
      
      ninew=max(round(ni2+(p1c-p1c2)*(ni2-ni1)/(p1c2-p1c1)),1)    
      ind=ind+1
      
      altdata=idata*(1+ninew*aa^2)
      
      
    bf0=BF(nulldata, n=ni1, m=mi, t=aa2) 
    bf1=BF(altdata, n=ni1, m=mi, t=aa2)
    
    k1=quantile(bf0, p0m)
    
    p1cnew=sum(bf1<(k1)) / length(bf1)
      
      print(paste0("It ",ind,": Suggested (n,m)=(", ninew,",",mi,") and achieved p1c=",p1cnew))
      
      vn=c(ni1,ni2,ninew) 
      if (length(unique(vn)) == 1) { 
        return(hnew)
      }
      
      while( abs(p1cnew-p1c)/p1c>0.001) {
        ni1=ni2
        p1c1=p1c2
        ni2=ninew
        p1c2=p1cnew
        
        while(p1c2==p1c1) {
          ni2=ni2+1
          ind=ind+1
          
          altdata=idata*(1+ni2*aa^2)
          
          
          bf0=BF(nulldata, n=ni1, m=mi, t=aa2)   
          bf1=BF(altdata, n=ni1, m=mi, t=aa2)
          
          k1=quantile(bf0, p0m)
    
          p1c2=sum(bf1<(k1)) / length(bf1)
          
          print(paste0("It ",ind,": Suggested (n,m)=(", ni2,",",mi,") and achieved p1c=",p1c2))
        }
        
        
        ninew=max(round(ni2+(p1c-p1c2)*(ni2-ni1)/(p1c2-p1c1)),1)    
        ind=ind+1
        
        altdata=idata*(1+ninew*aa^2)
        
        
        bf0=BF(nulldata, n=ni1, m=mi, t=aa2)   
        bf1=BF(altdata, n=ni1, m=mi, t=aa2)
        
        k1=quantile(bf0, p0m)
    
        p1cnew=sum(bf1<(k1)) / length(bf1)
        
        print(paste0("It ",ind,": Suggested (n,m)=(", ninew,",",mi,") and achieved p1c=",p1cnew))
        
        
        vn=c(ni1,ni2,ninew) 
        if (length(unique(vn)) < 3) { 
          break
        }
        
        
        
      }
    }
    
    print("End of Regula Falsi")
    #########################################################################################
    
    if (p1cnew < p1c) {
      
      while (p1cnew < p1c) {
        
        ninew=ninew+1
        ind=ind+1
        
        altdata=idata*(1+ninew*aa^2)
        
        
        bf0=BF(nulldata, n=ni1, m=mi, t=aa2)   
        bf1=BF(altdata, n=ni1, m=mi, t=aa2)
        
        k1=quantile(bf0, p0m)
    
        p1cnew=sum(bf1<(k1)) / length(bf1)
        
        print(paste0("It ",ind,": Suggested (n,m)=(", ninew,",",mi,") and achieved p1c=",p1cnew))
        
        
      } } else {
        
        while (p1cnew >= p1c) {
          
          niprev=ninew
          p1cprev=p1cnew
          
          ninew=ninew-1
          ind=ind+1
          
          altdata=idata*(1+ninew*aa^2)
          
          
          bf0=BF(nulldata, n=ni1, m=mi, t=aa2)  
          bf1=BF(altdata, n=ni1, m=mi, t=aa2)
          
          k1=quantile(bf0, p0m)
    
          p1cnew=sum(bf1<(k1)) / length(bf1)
          
          print(paste0("It ",ind,": Suggested (n,m)=(", ninew,",",mi,") and achieved p1c=",p1cnew))
          
        } 
        
        ninew=niprev
        p1cnew=p1cprev
        
      }
    
    RESULTS[ind0,1]=ninew
    RESULTS[ind0,3]=p1cnew
    RESULTS[ind0,4]=round(k1,3)

    print(RESULTS[1:ind0,])
    
  }
  
  return(RESULTS)
  
}


n_cond()









