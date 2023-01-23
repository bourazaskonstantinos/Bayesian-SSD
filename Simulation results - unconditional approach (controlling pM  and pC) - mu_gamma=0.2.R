
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

### a function that given pm and pc, returns the pairs of n and m  ### 


nm_pairs=function(pm, pc, mm = 3:17){
  
  RESULTS=data.frame(n=NA, m=mm, pc=NA, rk1=NA, k0=NA)
  
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
    ni1=round(500*pc^2/mi)
    
    # generate Q statistics under M_0
    
    set.seed(3)
    nulldata=rchisq(it, df=mi-1)
    
    # generate Q statistics under M_1
    set.seed(4)
    idata=rchisq(it, df=mi-1)
    altdata=idata*(1+ni1*aa^2)
    
    
    bf0=BF(nulldata, n=ni1, m=mi, t=aa2)
    bf1=BF(altdata, n=ni1, m=mi, t=aa2)
    
    k0=quantile(bf1, 1-pm)
    k1=quantile(bf0, pm)
    
    p1c=sum(bf1<(k1)) / length(bf1)
    p0c = sum(bf0 > k0) / length(bf0) # probability of decisive correct evidence under M_0

    pc1=(p0c+p1c)/2
    
    print(paste0("It ",ind,": Suggested (n,m)=(", ni1,",",mi,") and achieved pc=",pc1))
    
    if (abs(pc1-pc)/pc<0.001) {stop=1 ; ninew=ni1 ; pcnew=pc1}
    
    
    if (stop == 0)  {  
      ## second initial value (n)
      ind=2
      ni2=2*ni1
      
      # generate Q statistics under M_0
      
      altdata=idata*(1+ni2*aa^2)
      
      
      bf0=BF(nulldata, n=ni2, m=mi, t=aa2)
      bf1=BF(altdata, n=ni2, m=mi, t=aa2)
      
      k0=quantile(bf1, 1-pm)
      k1=quantile(bf0, pm)
      
      p1c =sum(bf1<(k1)) / length(bf1)
      p0c = sum(bf0 > k0) / length(bf0) # probability of decisive correct evidence under M_0
      
      pc2=(p0c+p1c)/2
      
      print(paste0("It ",ind,": Suggested (n,m)=(", ni2,",",mi,") and achieved pc=",pc2))
            
      if (abs(pc2-pc)/pc<0.001) {stop=1 ; ninew=ni2 ; pcnew=pc2}
      
      
      #########################################################################################
      
      while(pc2==pc1) {
        ni2=ni2+1
        ind=ind+1
        # generate Q statistics under M_0
        
        altdata=idata*(1+ni2*aa^2)
        
        
        bf0=BF(nulldata, n=ni2, m=mi, t=aa2)
        bf1=BF(altdata, n=ni2, m=mi, t=aa2)
        
        k0=quantile(bf1, 1-pm)
        k1=quantile(bf0, pm)
        
        p1c=sum(bf1<(k1)) / length(bf1)
        p0c = sum(bf0 > k0) / length(bf0) # probability of decisive correct evidence under M_0
        
        pc2=(p0c+pc1)/2
        
        print(paste0("It ",ind,": Suggested (n,m)=(", ni2,",",mi,") and achieved pc=",pc2))
      }
      
      if (abs(pc2-pc)/pc<0.001) {stop=1 ; ninew=ni2 ; pcnew=pc2}
      
      
    }
    #########################################################################################
    
    if (stop==0) {
      
      ninew=max(round(ni2+(pc-pc2)*(ni2-ni1)/(pc2-pc1)),1)    
      ind=ind+1
      # generate Q statistics under M_0
      
      altdata=idata*(1+ninew*aa^2)
      
      
      bf0=BF(nulldata, n=ninew, m=mi, t=aa2)
      bf1=BF(altdata, n=ninew, m=mi, t=aa2)
      
      k0=quantile(bf1, 1-pm)
      k1=quantile(bf0, pm)
      
      p1c=sum(bf1<(k1)) / length(bf1)
      p0c = sum(bf0 > k0) / length(bf0) # probability of decisive correct evidence under M_0
      
      pcnew=(p0c+p1c)/2
      
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
          
          k0=quantile(bf1, 1-pm)
          k1=quantile(bf0, pm)
          
          p1c=sum(bf1<(k1)) / length(bf1)
          p0c = sum(bf0 > k0) / length(bf0) # probability of decisive correct evidence under M_0
          
          pc2=(p0c+p1c)/2
          
          print(paste0("It ",ind,": Suggested (n,m)=(", ni2,",",mi,") and achieved pc=",pc2))
        }
        
        
        ninew=max(round(ni2+(pc-pc2)*(ni2-ni1)/(pc2-pc1)),1)    
        ind=ind+1

        altdata=idata*(1+ninew*aa^2)
        
        
        bf0=BF(nulldata, n=ninew, m=mi, t=aa2)
        bf1=BF(altdata, n=ninew, m=mi, t=aa2)
        
        k0=quantile(bf1, 1-pm)
        k1=quantile(bf0, pm)
        
        p1c=sum(bf1<(k1)) / length(bf1)
        p0c = sum(bf0 > k0) / length(bf0) # probability of decisive correct evidence under M_0
        
        pcnew=(p0c+p1c)/2
        
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
        # generate Q statistics under M_0
        
        altdata=idata*(1+ninew*aa^2)
        
        
        bf0=BF(nulldata, n=ninew, m=mi, t=aa2)
        bf1=BF(altdata, n=ninew, m=mi, t=aa2)
        
        k0=quantile(bf1, 1-pm)
        k1=quantile(bf0, pm)
        
        p1c=sum(bf1<(k1)) / length(bf1)
        p0c = sum(bf0 > k0) / length(bf0) # probability of decisive correct evidence under M_0
        
        pcnew=(p0c+p1c)/2
        
        print(paste0("It ",ind,": Suggested (n,m)=(", ninew,",",mi,") and achieved pc=",pcnew))
        
        
      } } else {
        
        while (pcnew >= pc) {
          
          niprev=ninew
          pcprev=pcnew
          
          ninew=ninew-1
          ind=ind+1
          # generate Q statistics under M_0
          
          altdata=idata*(1+ninew*aa^2)
          
          
          bf0=BF(nulldata, n=ninew, m=mi, t=aa2)
          bf1=BF(altdata, n=ninew, m=mi, t=aa2)
          
          k0=quantile(bf1, 1-pm)
          k1=quantile(bf0, pm)
          
          p1c=sum(bf1<(k1)) / length(bf1)
          p0c = sum(bf0 > k0) / length(bf0) # probability of decisive correct evidence under M_0
          
          pcnew=(p0c+p1c)/2
          
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

pairs2a=nm_pairs(pm=0.01, pc=0.8)
pairs2a
pairs2b=nm_pairs(pm=0.05, pc=0.8)
pairs2b


pairs2c=nm_pairs(pm=0.01, pc=0.9)
pairs2c
pairs2d=nm_pairs(pm=0.05, pc=0.9)
pairs2d

### results ###

### pm=0.01, pc=0.8 ###

pairs2a$n
pairs2a$m
pairs2a$rk1
pairs2a$k0

### pm=0.05, pc=0.8 ###

pairs2b$n
pairs2b$m
pairs2b$rk1
pairs2b$k0

### pm=0.01, pc=0.9 ###

pairs2c$n
pairs2c$m
pairs2c$rk1
pairs2c$k0

### pm=0.05, pc=0.9 ###

pairs2d$n
pairs2d$m
pairs2d$rk1
pairs2d$k0


pairs2e
pairs2e$n
pairs2e$m
pairs2e$rk1
pairs2e$k0
cf(pairs2e$n,pairs2e$m)

pairs2f
pairs2f$n
pairs2f$m
pairs2f$rk1
pairs2f$k0
cf(pairs2f$n,pairs2f$m)







