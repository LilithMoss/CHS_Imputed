file.Fun <- function(sig.snp,sig.chr){
  sig.file <- rep(NA,length(sig.snp))
  
  for(j in 1:length(sig.snp)){
    chr=sig.chr[j]
    snp=sig.snp[j]
    
    for(i in 0:9){
      t <- readRDS(paste0(search_path,"G.chr",chr,".000",i))
      a <- snp %in% t[,1]
      print(a);print(i)
      if(a==TRUE) break
    }
    
    if(a==FALSE){
      for(i in 10:99){
        t <- readRDS(paste0(search_path,"G.chr",chr,".00",i))
        a <- snp %in% t[,1]
        print(a);print(i)
        if(a==TRUE) break
      }
    } 
    
    if(a==FALSE){
      for(i in 100:350){
        t <- readRDS(paste0(search_path,"G.chr",chr,".0",i))
        a <- snp %in% t[,1]
        print(a);print(i)
        if(a==TRUE) break
      }
    } 
    
    if(nchar(i)==1){
      out<-paste0("000",i)
    } else if(nchar(i)==2){
      out<-paste0("00",i)
    } else 
      out <- paste0("0",i)
    
    sig.file[j] <- out
    a=FALSE
    
  } #j Loop
  
  return(sig.file)
  
} #Function