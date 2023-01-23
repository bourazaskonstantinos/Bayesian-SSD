
### iterations ###
it=5*10^4

### analysis prior ###
ss=10^4
set.seed(1)
a2=rt(ss, df=4)
s=1/7
aa2=abs(s*a2) # the sample from the analysis prior

# design prior
set.seed(2)
a=rt(it, df=4)
mu=0.2 ; s=1/55
aa=abs(s*a+mu) # the sample from the design prior


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################

### a function that given p0m and p1dc, returns the pairs of n and m  ### 


nm_pairs=function(p0m, p1c, mm = 3:17){
  
  RESULTS=data.frame(n=NA, m=mm, p1c=NA, rk1=NA)
  
  ### nested function that returns the BFs ###
  BF=function(Qm, n, m, t) {
    
    bfa=c()
    for (i in 1:it) {
      
      BFa=(n^((m-1)/2))*exp(-Qm[i]/2)/mean( (1/n+t^2)^((1-m)/2) * exp(-Qm[i]/2*(1+n*t^2)^(-1) ) )
      bfa=c(bfa,BFa)
      
    }  
    
    return(bfa)
    
  }
  
  
    
  ### modified regula falsi process to obtain the pairs of n, m
  
  ind0=0
  for (mi in mm) {
    
    ind0=ind0+1
    stop=0
    ## first initial value (n)
    ind=1
    ni1=round(500*p1c^2/mi)
    
    # generate Q statistics under M_0
    
    set.seed(3)
    nulldata=rchisq(it, df=mi-1)
    
    # generate Q statistics under M_1
    set.seed(4)
    idata=rchisq(it, df=mi-1)
    altdata=idata*(1+ni1*aa^2)
    
    
    bf0=BF(nulldata, n=ni1, m=mi, t=aa2)
    bf1=BF(altdata, n=ni1, m=mi, t=aa2)
    
    k1=quantile(bf0, p0m)
    
    p1c1=sum(bf1<(k1)) / length(bf1)

    
    print(paste0("It ",ind,": Suggested (n,m)=(", ni1,",",mi,") and achieved p1c=",p1c1))
    
    if (abs(p1c1-p1c)/p1c<0.001) {stop=1 ; ninew=ni1 ; p1cnew=p1c1}
    
    
    if (stop == 0)  {  
      ## second initial value (n)
      ind=2
      ni2=2*ni1
      
      altdata=idata*(1+ni2*aa^2)
      
      
    bf0=BF(nulldata, n=ni2, m=mi, t=aa2)  
    bf1=BF(altdata, n=ni2, m=mi, t=aa2)
    
    k1=quantile(bf0, p0m)
    
    p1c2=sum(bf1<(k1)) / length(bf1)



      
      print(paste0("It ",ind,": Suggested (n,m)=(", ni2,",",mi,") and achieved p1c=",p1c2))
            
      if (abs(p1c2-p1c)/p1c<0.001) {stop=1 ; ninew=ni2 ; pd1cnew=p1c2}
      
      
      #########################################################################################
      
      while(p1c2==p1c1) {
        ni2=ni2+1
        ind=ind+1
        # generate Q statistics under M_0
        
        altdata=idata*(1+ni2*aa^2)
        
        
    bf0=BF(nulldata, n=ni2, m=mi, t=aa2)  
    bf1=BF(altdata, n=ni2, m=mi, t=aa2)
    
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
      # generate Q statistics under M_0
      
      altdata=idata*(1+ninew*aa^2)
      
      
    bf0=BF(nulldata, n=ninew, m=mi, t=aa2) 
    bf1=BF(altdata, n=ninew, m=mi, t=aa2)
    
    k1=quantile(bf0, p0m)
    
    p1cnew=sum(bf1<(k1)) / length(bf1)
      
      print(paste0("It ",ind,": Suggested (n,m)=(", ninew,",",mi,") and achieved p1c=",p1cnew))
      
      vn=c(ni1,ni2,ninew) 
      if (length(unique(vn)) == 1) { 
        message("Warning message: the desired accuracy cannot be achieved.")
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
          # generate Q statistics under M_0
          
          altdata=idata*(1+ni2*aa^2)
          
          
          bf0=BF(nulldata, n=ni2, m=mi, t=aa2)   
          bf1=BF(altdata, n=ni2, m=mi, t=aa2)
          
          k1=quantile(bf0, p0m)
    
          p1c2=sum(bf1<(k1)) / length(bf1)
          
          print(paste0("It ",ind,": Suggested (n,m)=(", ni2,",",mi,") and achieved p1c=",p1c2))
        }
        
        
        ninew=max(round(ni2+(p1c-p1c2)*(ni2-ni1)/(p1c2-p1c1)),1)    
        ind=ind+1

        altdata=idata*(1+ninew*aa^2)
        
        
        bf0=BF(nulldata, n=ninew, m=mi, t=aa2)   
        bf1=BF(altdata, n=ninew, m=mi, t=aa2)
        
        k1=quantile(bf0, p0m)
    
        p1cnew=sum(bf1<(k1)) / length(bf1)
        
        print(paste0("It ",ind,": Suggested (n,m)=(", ninew,",",mi,") and achieved p1c=",p1cnew))
        
        
        vn=c(ni1,ni2,ninew) 
        if (length(unique(vn)) < 3) { 
          message("Warning message: the desired accuracy in Regula Falsi cannot be achieved.")
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
        # generate Q statistics under M_0
        
        altdata=idata*(1+ninew*aa^2)
        
        
        bf0=BF(nulldata, n=ninew, m=mi, t=aa2)   
        bf1=BF(altdata, n=ninew, m=mi, t=aa2)
        
        k1=quantile(bf0, p0m)
    
        p1cnew=sum(bf1<(k1)) / length(bf1)
        
        print(paste0("It ",ind,": Suggested (n,m)=(", ninew,",",mi,") and achieved p1c=",p1cnew))
        
        
      } } else {
        
        while (p1cnew >= p1c) {
          
          niprev=ninew
          p1cprev=p1cnew
          
          ninew=ninew-1
          ind=ind+1
          # generate Q statistics under M_0
          
          altdata=idata*(1+ninew*aa^2)
          
          
          bf0=BF(nulldata, n=ninew, m=mi, t=aa2)  
          bf1=BF(altdata, n=ninew, m=mi, t=aa2)
          
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

pairs2a=nm_pairs(p0m=0.01, p1c=0.8)
pairs2a
pairs2b=nm_pairs(p0m=0.05, p1c=0.8)
pairs2b


pairs2c=nm_pairs(p0m=0.01, p1c=0.9)
pairs2c
pairs2d=nm_pairs(p0m=0.05, p1c=0.9)
pairs2d

### results ###

### p0m=0.01, p1c=0.8 ###

pairs2a$n
pairs2a$m
pairs2a$rk1

### p0m=0.05, p1c=0.8 ###

pairs2b$n
pairs2b$m
pairs2b$rk1

### p0m=0.01, p1c=0.9 ###

pairs2c$n
pairs2c$m
pairs2c$rk1

### p0m=0.05, p1c=0.9 ###

pairs2d$n
pairs2d$m
pairs2d$rk1







