# setwd("./research/kareem/test_runs/R/4.1.2/") #for superior
# for simulation

#library
rm(list=ls())
library(svines)
library(parallel)
library(dplyr)
library(gcmr)
library(evmix)

# function
cl <- makeCluster(detectCores() - 1)
clusterEvalQ(cl, {
  rm(list=ls());
  library(svines);
  library(dplyr)
})
starttime <- Sys.time()
# functions
svinePred <- function(data, udata, model, var, tp){
  pred <- svinecop_sim(1, rep=200, model = model, past=udata[1:(tp-1), var])[1,,]
  if (length(var)==1){ 
    if (var=="x"){
      val<-mean(quantile(data[,"x"],pred))
    }else{val<-mean(quantile(data[,"y"],pred))} 
  }else{
    val<-c(mean(quantile(data[,"x"],pred["x",])),mean(quantile(data[,"y"],pred["y",])), mean(quantile(data[,"z"],pred["z",])) )
  }
  return(val)
}
svineGC <- function(sub, v1, v2){
  gc <- log(sum((v1-sub)^2)/sum((v2-sub)^2))
  return(gc)
}
# parameter
#size <- c(400)

T2 <- 20
T1 <- T2 * 0.5
model <- c("A1" )  # ,"A2","B1","B2","B3","B4")
N <- 19                                   # number of bootstrap samples
order <- 1
it=1#  <- 1:1000
#trace(svines:::to_unif, quote( u <- ifelse(u<=0, 0, ifelse(u>=1, 1, u)) ), at=3, print=FALSE)
# simulation
result <- matrix(0, nrow=length(model), ncol=2, dimnames=list(model, c("BETAyx","yz.x")))
for(j  in 1:length(model)){
  dgp <- model[j]
  #T2 <- size[j]
  #T1 <- as.integer(0.5*T2)
  pmat <- matrix(0, nrow=length(it), ncol=2, dimnames=list(it, c("BETAyx","yz.x")))
  
for(i in it){
    if(i %% 10 == 1) cat("Model", dgp, "T2=", T2, "it=", i, "\n")
    
    ###########################################################################################################
    # data generate
    if(dgp=="A1"){
      y <- arima.sim(list(ar=0.5), n=T2, rand.gen=rnorm)
      x <- arima.sim(list(ar=0.5), n=T2, rand.gen=rnorm)
      z <- arima.sim(list(ar=0.5), n=T2, rand.gen=rnorm)
    }
    if(dgp=="B1") {
      y <- arima.sim(list(ar=0.5), n=T2, rand.gen=rnorm)
      x <- rnorm(1)
      for(k  in 2:T2){
        x[k] <- abs(x[k-1])^0.8 + rnorm(1)
      }
      z <- arima.sim(list(ar=0.5), n=T2, rand.gen=rnorm)
    }
    data <- data.frame(x=x, y=y,z=z)
    udata <- pseudo_obs(data)
    #edit-- next two lines  
    # write.table(y, paste0("./finalY/",dgp,"/T2_",T2,"_it_",i,".txt"), row.names=FALSE, col.names=FALSE)
    #write.table(x, paste0("./finalX/",dgp,"/T2_",T2,"_it_",i,".txt"), row.names=FALSE, col.names=FALSE)
 #release write lines above; add z
    ###########################################################################################################
    # S-vine model
    Sx <- svinecop(udata[1:T1,'x'], p=order)
    S2 <- svinecop(udata[1:T1,], p=order)
    # predict
    predsX <- sapply(c(T1+1):T2, function(p) svinePred(data, udata, Sx, "x", p))
    preds2 <- sapply(c(T1+1):T2, function(p) svinePred(data, udata, S2, c("x","y","z"), p))
    rownames(preds2)<-c("x","y","z")
    # estimate GC statistics
    subs <- data[(T1+1):T2,] #(sub, v1, v2)
    GCyx_stat <- svineGC(subs[,'x'], predsX, preds2['x',])
    # generate the bootstrap samples from the null distribution
    # use M-vine to generate the bootstrap samples
    #cs_structure <- matrix(c(1,2,1,0), 2, 2)
    cs_structure <- matrix(c(2,1,3, 1,2,0, 1,0,0), 3, 3)
    #out_vertices <- in_vertices <- c(1,2)
    out_vertices <- in_vertices <- c(1,2,3)
    #S2_boot <- svinecop(udata[1:T1,], p=order, cs_structure=cs_structure, out_vertices=out_vertices, in_vertices=in_vertices)
    S2_boot <- svinecop(udata[1:T1,], p=order, cs_structure=cs_structure, out_vertices=out_vertices, in_vertices=in_vertices)
                        
    sampleX<-matrix(0,T2,N)
    sampleY<-matrix(0,T2,N)
    sampleZ<-matrix(0,T2,N)
    
    genYfromX<-S2_boot$pair_copulas[[1]][[1]] # generating Y_t-2 form X_t-2
    genXfromX<-S2_boot$pair_copulas[[1]][[2]] # generating X_t-1 form X_t-2 
    genZfromX<-S2_boot$pair_copulas[[1]][[3]] # why [[4]]
    for (iterN in 1:N){
      Umat<-matrix(runif(3*T2),T2,3)
      colnames(Umat)<-c("x","y","z")
      # current stop
      for (iter in 1:T2){
        Umat[iter,"y"]<-inverse_rosenblatt(Umat[iter,c("x","y")], genYfromX, cores = 1)[2]
        Umat[iter,"z"]<-inverse_rosenblatt(Umat[iter,c("x","z")], genZfromX, cores = 1)[2]
        if (iter<T2){
          Umat[iter+1,"x"]<-inverse_rosenblatt(matrix(Umat[iter:(iter+1),1],1,2), genXfromX, cores = 1)[2]}
      }
      sampleX[,iterN]<-as.vector(quantile(data[,"x"],Umat[,"x"]))
      sampleY[,iterN]<-as.vector(quantile(data[,"y"],Umat[,"y"]))
      sampleZ[,iterN]<-as.vector(quantile(data[,"z"],Umat[,"z"]))
    }
    # Here is what takes the most time to run ===========================      
    clusterEvalQ(cl, rm(list=ls()))
    clusterExport(cl, c("T1","T2","order","sampleX","sampleY","sampleZ","svinePred", "svineGC"))
    
    res_parLapply <- parLapply(cl, 1:N, function(n) {
      
      tryCatch(expr={
        
        #trace(svines:::to_unif, quote( u <- ifelse(u<=0, 0, ifelse(u>=1, 1, u)) ), at=3, print=FALSE)
        
        # resample series
        sample <- data.frame(x=sampleX[,n],y=sampleY[,n],z=sampleZ[,n]) 
        usample <- pseudo_obs(sample)
        
        # statistic
        Cx <- svinecop(usample[1:T1,'x'], p=order)
        C2 <- svinecop(usample[1:T1,], p=order)
        
        # predict
        predX <- sapply(c(T1+1):T2, function(p) svinePred(sample, usample, Cx, "x", p))
        pred2 <- sapply(c(T1+1):T2, function(p) svinePred(sample, usample, C2, c("x","y","z"), p))
        rownames(pred2)<-c("x","y","z")
        #untrace(svines:::to_unif)
        # estimate GC statistics
        sub <- sample[(T1+1):T2,]
        #GCyx <- svineGC(subs[,'x'], predsX, preds2['x',])
        GCyz.x <- svineGC(sub[,'x'], predX, pred2['x',])
        GCyz.x
      }, error=function(e) NA, warning=function(w) NA)
    })
    # to here ===============================    (runs for about 7 minutes )  
    #hist(unlist(res_parLapply),xlim=range(c(unlist(res_parLapply),GCyx_stat)))
    #points(GCyx_stat, 0, cex=1, pch=16,col=2)
    if(sum(is.na(res_parLapply)) >= N*0.1) break
    pmat[i, "yz.x"] <- mean(unlist(na.omit(res_parLapply)) >= GCyx_stat)
    #cat("it=", i,"reject or not=",as.numeric(pmat[i,]<0.05),"\n")
    write.csv(as.data.frame(pmat), paste0("./pvals_T2_",T2,"dgp",dgp,".csv"))
    
  } #closing for i loop
  result[j,] <- apply(pmat,2,function(x) mean(x<0.05))
write.csv(as.data.frame(result), "./result.csv")
write.csv(as.data.frame(result), paste0("./results",dgp,"_T2_",T2,".csv") )
} #closing for j loop

stopCluster(cl)
#untrace(svines:::to_unif)
endtime <- Sys.time()
starttime
endtime
duration <- (as.numeric(endtime) - as.numeric(starttime)) / 3600 # duration in hours
duration
result
