evd <- function(Q){  
  decompo <- eigen(Q)
  lambda <- decompo$values
  GAMMA <- decompo$vectors
  invGAMMA <- solve(GAMMA)
  list(GAMMA=GAMMA, invGAMMA=invGAMMA, lambda=lambda)
}  

# without the diag, quite a bit faster
#expm_evd <- function(x, el) x$GAMMA %*% diag(exp(x$lambda * el)) %*% x$invGAMMA
expm_evd <- function(x, el)x$GAMMA %*% (exp(x$lambda * el) * x$invGAMMA)


## function for conditional likelihoods at nodes
## written by Liam J. Revell 2015, 2016, 2019, 2020
## with input from (& structural similarity to) function ace by E. Paradis et al. 2013
## using eigen value decomposition for symmetric model
fitMk2 <- function(tree,x,model="SYM",fixedQ=NULL,...){
  if(hasArg(output.liks)) output.liks<-list(...)$output.liks
  else output.liks<-FALSE
  if(hasArg(q.init)) q.init<-list(...)$q.init
  else q.init<-length(unique(x))/sum(tree$edge.length)
  if(hasArg(opt.method)) opt.method<-list(...)$opt.method
  else opt.method<-"nlminb"
  if(hasArg(min.q)) min.q<-list(...)$min.q
  else min.q<-1e-12
  N<-Ntip(tree)
  M<-tree$Nnode
  if(model=="SYM" || model=="ER") is_symm <- TRUE 
  else is_symm <- FALSE  
  if(is.matrix(x)){
    x<-x[tree$tip.label,]
    m<-ncol(x)
    states<-colnames(x)
  } else {
    x<-to.matrix(x,sort(unique(x)))
    x<-x[tree$tip.label,]
    m<-ncol(x)
    states<-colnames(x)
  }
  if(hasArg(pi)) pi<-list(...)$pi
  else pi<-"equal"
  if(is.numeric(pi)) root.prior<-"given"
  if(pi[1]=="equal"){ 
    pi<-setNames(rep(1/m,m),states)
    root.prior<-"flat"
  } else if(pi[1]=="estimated"){ 
    pi<-if(!is.null(fixedQ)) statdist(fixedQ) else 
      statdist(summary(fitMk(tree,x,model),quiet=TRUE)$Q)
    cat(paste("Using pi estimated from the stationary",
              "distribution of Q assuming a flat prior.\npi =\n"))
    print(round(pi,6))
    cat("\n")
    root.prior<-"stationary"
  } else if(pi[1]=="fitzjohn") root.prior<-"nuisance"
  if(is.numeric(pi)){ 
    pi<-pi/sum(pi)
    if(is.null(names(pi))) pi<-setNames(pi,states)
    pi<-pi[states]
  } 
  if(is.null(fixedQ)){
    if(is.character(model)){
      rate<-matrix(NA,m,m)
      if(model=="ER"){ 
        k<-rate[]<-1
        diag(rate)<-NA
      } else if(model=="ARD"){
        k<-m*(m-1)
        rate[col(rate)!=row(rate)]<-1:k
      } else if(model=="SYM"){
        k<-m*(m-1)/2
        ii<-col(rate)<row(rate)
        rate[ii]<-1:k
        rate<-t(rate)
        rate[ii]<-1:k
      }
    } else {
      if(ncol(model)!=nrow(model)) 
        stop("model is not a square matrix")
      if(ncol(model)!=ncol(x)) 
        stop("model does not have the right number of columns")
      rate<-model
      k<-max(rate)
    }
# need Q      
    Q<-matrix(1,m,m)
    diag(Q)<-0
    diag(Q)<- -rowSums(Q)
  } else {
    rate<-matrix(NA,m,m)
    k<-m*(m-1)
    rate[col(rate)!=row(rate)]<-1:k
    Q<-fixedQ
  }
  index.matrix<-rate
  tmp<-cbind(1:m,1:m)
  rate[tmp]<-0
  rate[rate==0]<-k+1
  liks<-rbind(x,matrix(0,M,m,dimnames=list(1:M+N,states)))
  pw<-reorder(tree,"postorder")
  lik<-function(Q,output.liks=FALSE,pi,...){
    if(hasArg(output.pi)) output.pi<-list(...)$output.pi
    else output.pi<-FALSE
    if(is.Qmatrix(Q)) Q<-unclass(Q)
    if(any(is.nan(Q))||any(is.infinite(Q))) return(1e50)
    if(is_symm) Q_evd <- evd(Q)
    comp<-vector(length=N+M,mode="numeric")
    parents<- pw$edge[,1]
    desc<-pw$edge[,2]
    el<-pw$edge.length
    root <- min(parents)
    vv <- rep(1, m)
    anc <- parents[1]
    for(i in 1:length(parents)){
      if(anc==parents[i]){
        if(is_symm) vv <- vv * expm_evd(Q_evd,el[i])%*%liks[desc[i],]
        else vv <- vv * expm(Q*el[i])%*%liks[desc[i],]
      }
      else {
        comp[anc]<-sum(vv)
        liks[anc,]<-vv/comp[anc]
        if(is_symm) vv <- expm_evd(Q_evd,el[i])%*%liks[desc[i],]
        else vv <- expm(Q*el[i])%*%liks[desc[i],]
        anc <- parents[i]
      }
    }
    if(is.numeric(pi)) vv<- vv*pi
    else if(pi[1]=="fitzjohn"){
      pi<-vv/sum(vv)
      vv<-vv*vv/sum(vv)
    }
    comp[anc]<-sum(vv)
    liks[anc,]<-vv/comp[anc]
    
    if(output.liks) return(liks[1:M+N,,drop=FALSE])
    else if(output.pi) return(pi)
    else {
      logL<--sum(log(comp[1:M+N]))
      if(is.na(logL)) logL<-Inf
      return(logL)
    }
  }
  if(is.null(fixedQ)){
    if(length(q.init)!=k) q.init<-rep(q.init[1],k)
    if(opt.method=="optim")
      fit<-optim(q.init,function(p) lik(makeQ(m,p,index.matrix),pi=pi),
                 method="L-BFGS-B",lower=rep(min.q,k))
    else if(opt.method=="none")
      fit<-list(objective=lik(makeQ(m,q.init,index.matrix),pi=pi),
                par=q.init)
    else	
      fit<-nlminb(q.init,function(p) lik(makeQ(m,p,index.matrix),pi=pi),
                  lower=rep(0,k),upper=rep(1e50,k))
    if(pi[1]=="fitzjohn") pi<-setNames(
      lik(makeQ(m,fit$par,index.matrix),FALSE,pi=pi,output.pi=TRUE),
      states)
    obj<-list(logLik=
                if(opt.method=="optim") -fit$value else -fit$objective,
              rates=fit$par,
              index.matrix=index.matrix,
              states=states,
              pi=pi,
              method=opt.method,
              root.prior=root.prior)
    if(output.liks) obj$lik.anc<-lik(makeQ(m,obj$rates,index.matrix),TRUE,
                                     pi=pi)
  } else {
    fit<-lik(Q,pi=pi)
    if(pi[1]=="fitzjohn") pi<-setNames(lik(Q,FALSE,pi=pi,output.pi=TRUE),states)
    obj<-list(logLik=-fit,
              rates=Q[sapply(1:k,function(x,y) which(x==y),index.matrix)],
              index.matrix=index.matrix,
              states=states,
              pi=pi,
              root.prior=root.prior)
    if(output.liks) obj$lik.anc<-lik(makeQ(m,obj$rates,index.matrix),TRUE,
                                     pi=pi)
  }
  lik.f<-function(q) -lik(q,output.liks=FALSE,
                          pi=if(root.prior=="nuisance") "fitzjohn" else pi)
  obj$lik<-lik.f
  class(obj)<-"fitMk"
  return(obj)
} 


rerootingMethod2 <- function(tree,x,model=c("ER","SYM"),...){
  if(!inherits(tree,"phylo")) 
    stop("tree should be an object of class \"phylo\".")
  if(hasArg(tips)) tips<-list(...)$tips
  else tips<-NULL
  if(!is.matrix(model)) model<-model[1]
  n<-Ntip(tree)
  # if vector convert to binary matrix
  if(!is.matrix(x)){ 
    yy<-to.matrix(x,sort(unique(x)))
    if(is.null(tips)) tips<-FALSE
  } else { 
    if(is.null(tips)) tips<-TRUE
    yy<-x
  }
  yy<-yy[tree$tip.label,]
  yy<-yy/rowSums(yy)
  YY<-fitMk2(tree,yy,model=model,output.liks=TRUE,...)
  Q<-matrix(c(0,YY$rates)[YY$index.matrix+1],length(YY$states),
            length(YY$states),dimnames=list(YY$states,YY$states))
  diag(Q)<--colSums(Q,na.rm=TRUE)
  if(isSymmetric(Q)){
    Q_evd <- evd(Q)
    is_symm <- TRUE
  }
  else is_symm <- FALSE

  tree <- reorder(tree, "cladewise")
  
  parents <- tree$edge[,1]
  desc <- tree$edge[,2]
  el <- tree$edge.length
  
  XX <- rbind(yy, YY$lik.anc)
  
  for (i in seq_along(el)) {
    if (child[j] > nTips) {
      if(is_symm) P <- expm_evd(Q_evd, el[i])
      else P <- expm(Q*el[i])
      tmp <- XX[parent[i],] / (XX[desc[i],] %*% P)
      XX[desc[i],] <- (tmp %*% P) * XX[desc[i],]
    }
  }
  
  #  if(tips) XX<-rbind(XX[1:n,],YY$lik.anc[1,], if(tree$Nnode>1) 
  #    XX[(n+1):nrow(XX),])
  #  else XX<-rbind(yy,YY$lik.anc[1,],if(tree$Nnode>1) XX)
  #  rownames(XX)<-1:(tree$Nnode+n)
  #  if(tips) rownames(XX)[1:n]<-tree$tip.label
  #  XX<-if(tips) XX else XX[1:tree$Nnode+n,]
  XX<-if(!tips) XX[1:tree$Nnode+n,]
  obj<-list(loglik=YY$logLik,Q=Q,marginal.anc=XX,tree=tree,x=yy)
  class(obj)<-"rerootingMethod"
  obj
}
