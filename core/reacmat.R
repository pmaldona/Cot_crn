library("limSolve")
library("XML")
library("rlist")
Rcpp::sourceCpp('./core/rn.cpp')

rn.xmlextract <- function(f, verbose=F) {
  p <- xmlRoot(xmlParse(f))[["model"]]
  #p <<- p
  sp.id <- NULL
  sp.name <- NULL
  if ("listOfSpecies" %in% names(p)) {
    sp.l <- xmlApply(p[["listOfSpecies"]],xmlAttrs)
    sp.id <- as.vector(sapply(sp.l,function(e) e["id"]))
    sp.name <- as.vector(sapply(sp.l,function(e) e["name"]))
  }
  sp.idn <- sp.id
  i <- which(grepl("^S[0-9]+$",sp.idn))
  sp.idn[i] <- sp.name[i]
  p <- p[["listOfReactions"]]
  n <- which(xmlSApply(p, xmlName) == "reaction")
  a <- xmlApply(p, xmlAttrs)
  #a <- as.list(as.data.frame(xmlSApply(p, xmlAttrs)))
  l <- list()
  n <- which(xmlApply(p, xmlName) == "reaction")
  R <- list()
  # print(a)
  for (i in n) {
    reversible <- a[[i]]["reversible"]
    # print(a[[i]])
    # print(reversible)
    reversible <- if (!is.na(reversible) && reversible=="false") F else T
    #print(reversible)
    reactants <- if("listOfReactants" %in% names(p[[i]])) xmlApply(p[[i]][["listOfReactants"]], xmlAttrs) else NULL
    products <- if("listOfProducts" %in% names(p[[i]])) xmlApply(p[[i]][["listOfProducts"]], xmlAttrs) else NULL
    r.stoichiometry <- sapply(reactants,function(e) { s <- e["stoichiometry"]; if (is.na(s)) 1 else s } )
    p.stoichiometry <- sapply(products,function(e) { s <- e["stoichiometry"]; if (is.na(s)) 1 else s } )
    names(r.stoichiometry) <- NULL; names(p.stoichiometry) <- NULL
    reactants <- if(is.null(reactants)) NULL else sapply(reactants,function(e) e["species"])
    products <- if(is.null(products)) NULL else sapply(products,function(e) e["species"])
    names(reactants) <- NULL; names(products) <- NULL
    R <- c(R,list(list( reactants=reactants,
                        products=products,
                        r.stoichiometry=as.numeric(r.stoichiometry),
                        p.stoichiometry=as.numeric(p.stoichiometry) )))
    if (reversible)
      R <- c(R,list(list( reactants=products,
                          products=reactants,
                          r.stoichiometry=as.numeric(p.stoichiometry),
                          p.stoichiometry=as.numeric(r.stoichiometry) )))
    if (verbose) {
      print("----------------------------")
      print(reversible)
      print(c(reactants,"-->",products))
      print(c(r.stoichiometry,"-->",p.stoichiometry))
    }
  }
  list(sp.id=sp.id,sp.name=sp.name,sp.idn=sp.idn,reac=R)
}

rn.proc <- function(f=rn.test,tot_org=T) {
  if (is.numeric(f)) {
    f <- paste0("./ReacNet/",dir("ReacNet","*.xml")[f])
  } 
  cat("f =",f,"\n")
  rn <- rn.xmlextract(f)
  scod <- 1:length(rn$sp.id); names(scod) <- rn$sp.id
  mr <- mp <- matrix(0,length(scod),length(rn$reac))
  rownames(mr) <- rn$sp.id; rownames(mp) <- rn$sp.id

  
  i <- 1
  for (r in rn$reac) {
    l <- length(r$reactants); lp <- length(r$products)
    j <- 1; for (s in r$reactants) { mr[s,i] <- r$r.stoichiometry[j]; j <- j + 1 } 
    j <- 1; for (s in r$products) { mp[s,i] <- r$p.stoichiometry[j]; j <- j + 1 } 
    i <- i + 1
  }
  cat(nrow(mr),"especies,",ncol(mr), "reacciones\n")
  #rn <<- c(rn,list(spc=scod,mr=mr,mp=mp))
  
  if(!tot_org){
    rn <- c(rn,list(spc=scod,mr=mr,mp=mp))
    nsp <- rn.linp_org(rn,rn$sp.id,F)
    rn$spc=scod
    rn$mr=mr
    rn$mp=mp
    rn$nsp=nsp
    rn <<- rn
    return(rn)
  }
  # return(rn)
  m=mp-mr
  
  rsp <- which(sapply(1:length(m[,1]),function(i) any(m[i,]!=0))) # especies reactivas
  csp <- which(sapply(1:length(m[,1]),function(i) all(m[i,]<=0))) # especies consumidas o no reactivas
  rcsp <- intersect(rsp,csp) # especies consumidas y reactivas
  if(length(rcsp!=0)){
    cm <- matrix(0,length(scod),length(rcsp))
    mr <- cbind(mr,cm)
    
    for(i in 1:length(rcsp)){
      cm[rcsp[i],i]=1
      rn$reac <- list.append(rn$reac,list(reactants=NULL,products=rcsp[i],r.stoichiometry=numeric(0),p.stoichiometry=1))
    } 
    mp <- cbind(mp,cm)
    
    
  }
  
  
  rn$spc=scod
  rn$mr=mr
  rn$mp=mp
  nsp <- rn.linp_org(rn,rn$sp.id,F)
  
  if(length(nsp)==0) {
    rn <<- rn
    return(rn)
  }
  
  cm <- matrix(0,length(scod),length(nsp))
  rn$mr <- cbind(rn$mr,cm)
  for(i in 1:length(nsp)){ 
    cm[nsp[i],i]=1
    rn$reac <- list.append(rn$reac,list(reactants=NULL,products=nsp[i],r.stoichiometry=numeric(0),p.stoichiometry=1))
  }
  rn$mp <- cbind(rn$mp,cm)
  
  rn <<- rn
  return(rn)
}

rn.ccod <- function(rn) { 
  
  scod <- 1:length(rn$sp.id); names(scod) <- rn$sp.id
  rn <- rn$reac
  
  fn <- function(r){
    c(list(reactants = 0+scod[r$reactants], products = 0+scod[r$products]), r[3:4])
    
  }
  crn <- lapply(rn,fn)
  rn_init(crn,length(scod),F)
}

rn.display <- function(rn,i=1:ncol(rn$mr)) {
  for (k in i) {
    k.r <- which(rn$mr[,k]>0)
    k.p <- which(rn$mp[,k]>0)
    m <- cbind(rbind(rn$mr[k.r,k],rn$sp.idn[k.r]),rbind(rn$mp[k.p,k],rn$sp.idn[k.p]))
    m <- rbind(m," + ")
    if(length(k.r)>0) m[3,length(k.r)] <- " --> "
    else m <- cbind(c("","∅"," --> "),m)
    if(length(k.p)==0) m <- cbind(m,c("","∅",""))
    m[3,ncol(m)] <- "\n"
    m[1,m[1,]=="1"] <- ""
    cat(k,"\t",m,sep="")
  }
}

rn.mdraw <- function(m) {
  col <- colorRampPalette(c("darkblue","lightblue"))(4)
  col <- c(col,colorRampPalette(c("lightblue1","white","lightpink1"))(3))
  col <- c(col,colorRampPalette(c("red","yellow"))(4))
  breaks <- c(c(-1000,-4:-1)+.01,-.01,.01,c(1:4,1000)-.01)
  image(1:ncol(m),1:nrow(m),t(m),xlim=c(.5,ncol(m)+.5),ylim=c(.5,nrow(m)+.5),col=col,breaks=breaks,useRaster=T)
}

rn.linp_org <- function(crn,id=crn$sp.id,verbose=T) {
  Ns <- nrow(crn$mr); Nr <- ncol(crn$mr)
  E <- cbind(diag(Ns),-diag(Ns),crn$mp-crn$mr) # se agregan reacciones virtuales de introducción y destrucción de especies
  f <- rep(0,Ns) # la producción se iguala a 0 (dado que agregamos reacciones virtuales de destrucción)
  G <- diag(ncol(E))
  h <- rep(1.0,ncol(E)); h[1:(2*Ns)] <- 0 # las reacciones virtuales pueden tener tasa 0, las otras mayores que 0 (arbitrariamente 1.0)
  Cost <- rep(0,ncol(E)); Cost[1:Ns] <- 1 # idealmente costo 0 (ninguna reacción virtual de introducción operando)
  rn.linp.r <<- linp(E,f,G,h,Cost)
  X <- rn.linp.r$X
  k <- which(X[(1:Ns)]>0)
  out <- k
  if (verbose) {
    cat("especies necesarias: ",id[k],"\n")
    k <- which(X[(1:Ns)+Ns]>0)
    if (length(k)>0) cat("especies sobreproducidas adicionales: ",id[k],"\n")
  }
  return(out)
}

rn.linp_sp <- function(crn,id=crn$sp.id,sp.add.i=integer(0),verbose=F) {
  Ns <- nrow(crn$mr); Nr <- ncol(crn$mr)
  E <- cbind(diag(Ns),crn$mp-crn$mr) # se agregan reacciones virtuales de introducción de especies
  f <- rep(1,Ns) # la producción se iguala a 1 (siempre posible dado que agregamos reacciones virtuales de introducción)
  Cost <- rep(0,ncol(E)); Cost[setdiff(1:Ns,sp.add.i)] <- 1 # costo minimiza reacciones virtuales de introducción operando
  rn.linp.r <- linp(E,f,Cost=Cost)
  rn.linp.r <<- rn.linp.r
  X <- rn.linp.r$X
  r <- union(which(X[1:Ns]<1),sp.add.i)
  if (verbose) {
    cat("especies sobreproducibles: ",id[r],"\n")
  }
  r
}

rn.linp_org_ov <- function(crn,id=crn$sp.id,opsp,verbose=T) {
  Ns <- nrow(crn$mr); Nr <- ncol(crn$mr)
  
  O <- matrix(0,nrow=length(id),ncol=length(opsp))
  k=1
  for(i in opsp){ 
    O[i,k] <- -1
    k <- k+1
  }
  E <- cbind(O,crn$mp-crn$mr) # se agregan reacciones virtuales de introducción y destrucción de especies
  f <- rep(0,nrow(E)) # la producción se iguala a 0 (dado que agregamos reacciones virtuales de destrucción)
  G <- diag(ncol(E))
  h <- rep(1.0,ncol(E)) # las reacciones virtuales pueden tener tasa 0, las otras mayores que 0 (arbitrariamente 1.0)
  Cost <- rep(1.0,ncol(E)); # idealmente costo 0 (ninguna reacción virtual de introducción operando)
  rn.linp.r <<- linp(E,f,G,h,Cost)
  
  X <- rn.linp.r$X[-(1:length(opsp))]
  if(rn.linp.r$IsError) X <- NULL
  
  return(X)
}

rn.explore <- function(i=1:24,sp=T,sp.add=F) {
  for (k in i) {
    cat("--------------------------------\n",k,": ",sep="")
    rn.proc(k)
    rn.linp_org(rn)
    if (sp) {
      if (sp.add) sp.i <- which(rn.linp.r$X[(1:length(rn$sp.id))]>0) else sp.i <- integer(0)
      sp.op <- rn.sprod(rn,sp.i=sp.i,verbose=T)
      #sp.op <- rn.linp_sp(rn,crn$sp.id,sp.i)
      cat("especies sobreproducibles: ",rn$sp.id[sort(sp.op)],"\n")
      sp.nop <<- setdiff(1:length(rn$sp.id),sp.op)
      cat("especies no sobreproducibles: ",rn$sp.id[sp.nop],"\n")
      sp.op <<- sp.op
      sp.i <<- sp.i
    }
  }
}

rn.sprod <- function(crn,i=1:nrow(crn$mr),sp.i=integer(0),verbose=F) {
  Ns <- length(i); Nr <- ncol(crn$mr[i,])
  rc <- c()
  for(j in 1:ncol(crn$mr)) if(all(crn$mp[-i,j]==0) && all(crn$mr[-i,j]==0)) rc <- c(rc,j) 
  sp <- rep(F,Ns)
  if (length(sp.i)>0) {
    Esp <- numeric(Ns); Esp[sp.i] <- 1; # las especies declaradas a priori como sopreproducidas (sp.i) se introducen
    sp[sp.i] <- T # solo se analizarán las especies no declaradas sobreproducidas a priori
  }
  else Esp <- NULL
  E0 <- cbind(-diag(Ns),Esp,crn$mp[i,rc]-crn$mr[i,rc]) # se destruyen todas las especies (excepto la estudiada que se introducirá)
  f0 <- rep(0,Ns) # la producción se iguala a 0 (excepto para la especie estudiada que se igualará a 1)
  Cost0 <- rep(0,ncol(E0)) # costo 0, excepto para la especie estudiada que se introducirá con costo positivo
  for (s in which(!sp)) {
    if (verbose) cat("especie considerada:",crn$sp.id[s],"...")
    E <- E0; E[s,s] <- 1
    f <- f0; f[s] <- 1
    Cost <- Cost0; Cost[s] <- 1
    r <- linp(E,f,Cost=Cost)
    #if (s==36) { print(crn$sp.id[s]); E <<- E; r <<- r }
    if (r$solutionNorm==0) sp[s] <- T
    if (verbose) cat("costo:",r$solutionNorm,"\n")
    k <- which(r$X[1:Ns]>0)
    k <- k[k!=s]
    if (verbose && length(k)>0) cat("        sobreproducida:",crn$sp.id[k],"\n")
  }
  rn.linp.r <<- r
  if(any(sp))  return(i[which(sp)])
  else return(numeric(0))
}

rn.linp.flow <- function(crn) {
  Ns <- nrow(crn$mr); Nr <- ncol(crn$mr)
  E <- cbind(diag(Ns),-diag(Ns),crn$mp-crn$mr) # se agregan reacciones virtuales de introducción y destrucción de especies
  f <- rep(0,Ns) # la producción se iguala a 0 (dado que agregamos reacciones virtuales de destrucción)
  G <- diag(ncol(E))
  h <- rep(1.0,ncol(E)); h[1:(2*Ns)] <- 0 # las reacciones virtuales pueden tener tasa 0, las otras mayores que 0 (arbitrariamente 1.0)
  Cost <- rep(0,ncol(E)); Cost[1:Ns] <- 1 # idealmente costo 0 (ninguna reacción virtual de introducción operando)
  rn.linp.r <<- linp(E,f,G,h,Cost)
  if (verbose) {
    X <- rn.linp.r$X
    k <- which(X[(1:Ns)]>0)
    cat("especies necesarias: ",id[k],"\n")
    k <- which(X[(1:Ns)+Ns]>0)
    if (length(k)>0) cat("especies sobreproducidas adicionales: ",id[k],"\n")
  }
}

rn.dcom <- function(crn,sp){
  opsp <- rn.sprod(crn,sp)
  
  c_m <- crn$mr[sp,]
  nc_m <- c_m
  for(i in 1:ncol(crn$mr)){ 
    c_m[,i] <- (crn$mp[sp,i]!=0 & crn$mr[sp,i]!=0) & (crn$mp[sp,i]==crn$mr[sp,i])
    nc_m[,i] <-(crn$mp[sp,i]!=crn$mr[sp,i])
  }
  csp <- sp[rowSums(nc_m)==0 & rowSums(c_m)!=0]
  
  fsp <- setdiff(sp,union(opsp,csp))
  rc <- c()
  for(i in 1:ncol(crn$mr)) if(all(crn$mp[-sp,i]==0) && all(crn$mr[-sp,i]==0)) rc <- c(rc,i) 
  if(length(rc)==0) return()
  m <- crn$mp[fsp,rc,drop=F]-crn$mr[fsp,rc,drop=F]
  adj <- matrix(0,length(fsp),length(fsp))
  if(length(fsp)>0) for(i in 1:ncol(m)){
    # k.p <- which(m[,i]>0)
    # k.n <- which(m[,i]<0)
    k.a <- which(m[,i]!=0)
    adj[k.a,k.a] <- 1
    #if(length(k.p)>0 && length(k.n)>0) for(k in k.n) adj[k,k.p] <- 1 
  }
  adj <- adj+t(adj)
  adj[col(adj)==row(adj)] <-0 
  adj <- adj>0
  v <- fsp*0
  eci <- 0
  while(!is.na(ec <- match(0,v))){
    j <- 1
    while(T){
      ec<-c(ec,setdiff(which(adj[,ec[j]]),ec))
      j <- j+1
      if(j > length(ec)) break
    }
    eci <- eci+1
    v[ec] <- eci
  }
  return(list(opsp=opsp,csp=csp,fsp=fsp,ec=v))
}

rn.dcom.legacy <- function(crn,sp){
  opsp <- rn.sprod(crn,sp)
  fsp <- setdiff(sp,union(opsp,csp))
  rc <- c()
  for(i in 1:ncol(crn$mr)) if(all(crn$mp[-sp,i]==0) && all(crn$mr[-sp,i]==0)) rc <- c(rc,i) 
  if(length(rc)==0) return()
  m <- crn$mp[fsp,rc,drop=F]-crn$mr[fsp,rc,drop=F]
  adj <- matrix(0,length(fsp),length(fsp))
  if(length(fsp)>0) for(i in 1:ncol(m)){
    k.p <- which(m[,i]>0)
    k.n <- which(m[,i]<0)
    k.a <- which(m[,i]!=0)
    if(length(k.p)>0 && length(k.n)>0) for(k in k.n) adj[k,k.p] <- 1 
  }
  adj <- adj+t(adj)
  adj <- adj>0
  v <- fsp*0
  eci <- 0
  while(!is.na(ec <- match(0,v))){
    j <- 1
    while(T){
      ec<-c(ec,setdiff(which(adj[,ec[j]]),ec))
      j <- j+1
      if(j > length(ec)) break
    }
    eci <- eci+1
    v[ec] <- eci
  }
  return(list(opsp=opsp,fsp=fsp,ec=v))
}

rn.dcom.opsp <- function(crn,sp,opsp){
  
  c_m <- crn$mr[sp,]
  nc_m <- c_m
  for(i in 1:ncol(crn$mr)){ 
    c_m[,i] <- (crn$mp[sp,i]!=0 & crn$mr[sp,i]!=0) & (crn$mp[sp,i]==crn$mr[sp,i])
    nc_m[,i] <-(crn$mp[sp,i]!=crn$mr[sp,i])
  }
  csp <- sp[rowSums(nc_m)==0 & rowSums(c_m)!=0]
  
  fsp <- setdiff(sp,union(opsp,csp))
  rc <- c()
  
  for(i in 1:ncol(crn$mr)) if(all(crn$mp[-sp,i]==0) && all(crn$mr[-sp,i]==0)) rc <- c(rc,i) 
  if(length(rc)==0) return()
  m <- crn$mp[fsp,rc,drop=F]-crn$mr[fsp,rc,drop=F]
  adj <- matrix(0,length(fsp),length(fsp))
  if(length(fsp)>0) for(i in 1:ncol(m)){
    # k.p <- which(m[,i]>0)
    # k.n <- which(m[,i]<0)
    k.a <- which(m[,i]!=0)
    adj[k.a,k.a] <- 1
    #if(length(k.p)>0 && length(k.n)>0) for(k in k.n) adj[k,k.p] <- 1 
  }
  adj <- adj+t(adj)
  adj[col(adj)==row(adj)] <-0 
  adj <- adj>0
  v <- fsp*0
  eci <- 0
  while(!is.na(ec <- match(0,v))){
    j <- 1
    while(T){
      ec<-c(ec,setdiff(which(adj[,ec[j]]),ec))
      j <- j+1
      if(j > length(ec)) break
    }
    eci <- eci+1
    v[ec] <- eci
  }
  return(list(opsp=opsp,csp=csp,fsp=fsp,ec=v))
}

rn.dcom2sp <- function(rn,decom){
  opsp <- rn$sp.id[decom$opsp]
  csp <- rn$sp.id[decom$csp]
  fcs <- list()
  for(i in unique(decom$ec))
  {
    fcsp <- decom$fsp[decom$ec==i]
    # if (length(fcsp)==1) csp <- c(csp,rn$sp.id[fcsp])
    # else fcs <- c(fcs,list(rn$sp.id[fcsp]))
    fcs <- c(fcs,list(rn$sp.id[fcsp]))
  }
  return(list(opsp=opsp,fcs=fcs,csp=csp))
} 

# rn.adj <- function(crn,sp) {
#   m <- crn$mp-crn$mr
#   nsp <- nrow(m)
#   m[sp,] <- 0
#   m <- (m!=0)
#   clos <- function(s) {
#     b <- rowSums(m[,colSums(m[s,,drop=F])>0,drop=F])>0
#     s1 <- setdiff(which(b),s)
#     s <- union(s1,s)
#     while (length(s1)>0) {
#       b <- rowSums(m[,colSums(m[s1,,drop=F])>0,drop=F])>0
#       s1 <- setdiff(which(b),s)
#       s <- union(s1,s)
#     }
#     s
#   }
#   nz <- which(rowSums(m)!=0)
#   ceq <- numeric(nsp)
#   neq <- 1
#   for (i in nz) {
#     if (ceq[i]==0) {
#       ceq[clos(i)] <- neq
#       neq <- neq + 1
#     }
#   }
#   ceq
# }
# 
# rn.cfrag <- function(crn,sp) {
#   ceq <<- rn.adj(crn,sp)
#   N <- max(ceq)
#   cat(N,"circuitos frágiles detectados\n")
#   m <- crn$mr != 0
#   m[sp,] <- F
#   csm <- colSums(m)
#   if (N>0) for (i in 1:N) {
#     cat("circuito",i,"\n")
#     k <- which(ceq==i)
#     cat("especies:",crn$sp.idn[k],"\n")
#     cat("reacciones:\n")
#     j <- which(csm>0 && colSums(m[k,,drop=F])==csm)
#     rn.display(crn,j)
#   }
# }

rn.linp_rorg_c <- function(crn=crn_r,clos_idx,verbose=T) {
  r_m <- sapply(crn$s_reac,function(v)v)
  p_m <- sapply(crn$s_prod,function(v)v)
  sp <- crn$clos_sp[[clos_idx]]
  reac <- which(sapply(1:length(crn$reac),function(r) all(crn$reac[[r]] %in% sp)))
  mr <- r_m[sp,reac,drop=F]
  mp <- p_m[sp,reac,drop=F]
  if (verbose) { print(sp); print(reac); print(mr); print(mp) }
  Ns <- length(sp); Nr <- length(reac)
  E <- cbind(diag(Ns),-diag(Ns),mp-mr) # se agregan reacciones virtuales de introducción y destrucción de especies
  f <- rep(0,Ns) # la producción se iguala a 0 (dado que agregamos reacciones virtuales de destrucción)
  G <- diag(ncol(E))
  h <- rep(1.0,ncol(E)); h[1:(2*Ns)] <- 0 # las reacciones virtuales pueden tener tasa 0, las otras mayores que 0 (arbitrariamente 1.0)
  Cost <- rep(0,ncol(E)); Cost[1:Ns] <- 1 # idealmente costo 0 (ninguna reacción virtual de introducción operando)
  rn.linp.r <<- linp(E,f,G,h,Cost)
  if (verbose) {
    X <- rn.linp.r$X
    k <- which(X[(1:Ns)]>0)
    cat("especies necesarias: ",sp[k],"\n")
    k <- which(X[(1:Ns)+Ns]>0)
    if (length(k)>0) cat("especies sobreproducidas adicionales: ",sp[k],"\n")
  }
  return(all(rn.linp.r$X[1:Ns]==0))
}
# crn_r$ssm[(sapply(crn_r$ssm,function(i) rn.linp_org_c(crn_r,i,F)))]

rn.linp_corg_c <- function(crn=crn_r,clos_idx,verbose=T) {
  r_m <- sapply(crn$s_reac,function(v)v)
  p_m <- sapply(crn$s_prod,function(v)v)
  sp <- crn$ssmc_sp[[clos_idx]]
  reac <- which(sapply(1:length(crn$reac),function(r) all(crn$reac[[r]] %in% sp)))
  mr <- r_m[sp,reac,drop=F]
  mp <- p_m[sp,reac,drop=F]
  if (verbose) { print(sp); print(reac); print(mr); print(mp) }
  Ns <- length(sp); Nr <- length(reac)
  E <- cbind(diag(Ns),-diag(Ns),mp-mr) # se agregan reacciones virtuales de introducción y destrucción de especies
  f <- rep(0,Ns) # la producción se iguala a 0 (dado que agregamos reacciones virtuales de destrucción)
  G <- diag(ncol(E))
  h <- rep(1.0,ncol(E)); h[1:(2*Ns)] <- 0 # las reacciones virtuales pueden tener tasa 0, las otras mayores que 0 (arbitrariamente 1.0)
  Cost <- rep(0,ncol(E)); Cost[1:Ns] <- 1 # idealmente costo 0 (ninguna reacción virtual de introducción operando)
  rn.linp.r <<- linp(E,f,G,h,Cost)
  if (verbose) {
    X <- rn.linp.r$X
    k <- which(X[(1:Ns)]>0)
    cat("especies necesarias: ",sp[k],"\n")
    k <- which(X[(1:Ns)+Ns]>0)
    if (length(k)>0) cat("especies sobreproducidas adicionales: ",sp[k],"\n")
  }
  return(all(rn.linp.r$X[1:Ns]==0))
}

rn.linp_csorg_c <- function(crn=crn_r,clos_idx,verbose=T) {
  r_m <- sapply(crn$s_reac,function(v)v)
  p_m <- sapply(crn$s_prod,function(v)v)
  sp <- crn$clos_sp[[clos_idx]]
  reac <- which(sapply(1:length(crn$reac),function(r) all(crn$reac[[r]] %in% sp)))
  mr <- r_m[sp,reac,drop=F]
  mp <- p_m[sp,reac,drop=F]
  if (verbose) { print(sp); print(reac); print(mr); print(mp) }
  Ns <- length(sp); Nr <- length(reac)
  E <- cbind(diag(Ns),-diag(Ns),mp-mr) # se agregan reacciones virtuales de introducción y destrucción de especies
  f <- rep(0,Ns) # la producción se iguala a 0 (dado que agregamos reacciones virtuales de destrucción)
  G <- diag(ncol(E))
  h <- rep(1.0,ncol(E)); h[1:(2*Ns)] <- 0 # las reacciones virtuales pueden tener tasa 0, las otras mayores que 0 (arbitrariamente 1.0)
  Cost <- rep(0,ncol(E)); Cost[1:Ns] <- 1 # idealmente costo 0 (ninguna reacción virtual de introducción operando)
  rn.linp.r <<- linp(E,f,G,h,Cost)
  if (verbose) {
    X <- rn.linp.r$X
    k <- which(X[(1:Ns)]>0)
    cat("especies necesarias: ",sp[k],"\n")
    k <- which(X[(1:Ns)+Ns]>0)
    if (length(k)>0) cat("especies sobreproducidas adicionales: ",sp[k],"\n")
  }
  return(all(rn.linp.r$X[1:Ns]==0))
}

rn.n_reac_sp <- function(rn=crn_r) {
  k <- NULL
  for (l in rn$mgen) for (g in l) if (length(g)==1) k <- c(k,g)
  i <- 1:length(rn$s_reac[[1]])
  if (!is.null(k)) i <- i[-k]
  i
}

rn.void_inflow <- function(rn=crn_r) {
  for (l in rn$mgen) for (g in l) if (length(g)==0) return(0)
  return(1+length(rn.n_reac_sp(rn)))
}

rn.linp_sorg_c <- function(crn=crn_r,clos_idx,verbose=T) {
  r_m <- sapply(crn$s_reac,function(v)v)
  p_m <- sapply(crn$s_prod,function(v)v)
  sp <- crn$ssmo_sp[[clos_idx]]
  
  mr <- r_m[sp,reac,drop=F]
  mp <- p_m[sp,reac,drop=F]
  if (verbose) { print(sp); print(reac); print(mr); print(mp) }
  Ns <- length(sp); Nr <- length(reac)
  E <- cbind(diag(Ns),-diag(Ns),mp-mr) # se agregan reacciones virtuales de introducción y destrucción de especies
  f <- rep(0,Ns) # la producción se iguala a 0 (dado que agregamos reacciones virtuales de destrucción)
  G <- diag(ncol(E))
  h <- rep(1.0,ncol(E)); h[1:(2*Ns)] <- 0 # las reacciones virtuales pueden tener tasa 0, las otras mayores que 0 (arbitrariamente 1.0)
  Cost <- rep(0,ncol(E)); Cost[1:Ns] <- 1 # idealmente costo 0 (ninguna reacción virtual de introducción operando)
  rn.linp.r <<- linp(E,f,G,h,Cost)
  if (verbose) {
    X <- rn.linp.r$X
    k <- which(X[(1:Ns)]>0)
    cat("especies necesarias: ",sp[k],"\n")
    k <- which(X[(1:Ns)+Ns]>0)
    if (length(k)>0) cat("especies sobreproducidas adicionales: ",sp[k],"\n")
  }
  return(all(rn.linp.r$X[1:Ns]==0))
}

rn.n_reac_sp <- function(rn=crn_r) {
  k <- NULL
  for (l in rn$mgen) for (g in l) if (length(g)==1) k <- c(k,g)
  i <- 1:length(rn$s_reac[[1]])
  if (!is.null(k)) i <- i[-k]
  i
}

rn.void_inflow <- function(rn=crn_r) {
  for (l in rn$mgen) for (g in l) if (length(g)==0) return(0)
  return(1+length(rn.n_reac_sp(rn)))
}

# sum(sapply(crn_r$ssm,function(i) rn.linp_rorg_c(crn_r,i,F))) + rn.void_inflow()
# sum(sapply(1:length(crn_r$ssmc_sp),function(i) rn.linp_corg_c(crn_r,i,F))) + rn.void_inflow()

rn.is_elmnt <- function(l){
  if(length(l)==0) return(0)
  if(length(l)==1) return(1)
  f <- function(i){
    l_i <- l[-i]
    o <- l[[i]]
    k=which(sapply(l_i, function(b) all(is.element(b,o))))
    sapply(k, function(j) o <<- setdiff(o,l_i[[j]]))
    length(o)>0
  }
  sapply(1:length(l),f)
}

# discrete dynamic concentration calculation function
rn.dina_step <- function(crn,sp,x,v,n){
  
  # sp <- sp[x>0]
  rc <- c()
  for(i in 1:ncol(crn$mr)) if(all(crn$mp[-sp,i]==0) && all(crn$mr[-sp,i]==0)) rc <- c(rc,i) 
  C_M <- matrix(0,ncol=(n+1),nrow=length(x))
  rownames(C_M) <- crn$sp.id[sp]
  if(length(rc)==0){
    sapply(1:(n+1),function(i) C_M[,i] <- x)
    return(C_M)
  }
  S <- crn$mp[sp,rc,drop=F]-crn$mr[sp,rc,drop=F]
  
  Sm<- function(x){
    asp <- which(x>0)
    rc <- c()
    for(i in 1:ncol(crn$mr)) if(all(S[-asp,i]==0) && all(crn$mr[-asp,i]==0)) rc <- c(rc,i)
    St <- S
    St[,setdiff(1:ncol(S),rc)] <-0
    return(St)
  }
  for (i in 1:(n+1)){
    if(i==1) C_M[,i] <- x
    else C_M[,i] <- Sm(C_M[,i-1])%*%v + C_M[,i-1]
  }
  return(C_M)
  
}

# discrete dynamic concentration calculation function
rn.dina_step_stoi <- function(crn,sp,x,v,n){
  
  # sp <- sp[x>0]
  rc <- c()
  for(i in 1:ncol(crn$mr)) if(all(crn$mp[-sp,i]==0) && all(crn$mr[-sp,i]==0)) rc <- c(rc,i) 
  P_M <- matrix(0,ncol=(n),nrow=length(x))
  C_M <- matrix(0,ncol=(n+1),nrow=length(x))
  V_M <- matrix(0,ncol=(n+1),nrow=length(v))
  R_L <- list()
  rownames(C_M) <- crn$sp.id[sp]
  rownames(P_M) <- crn$sp.id[sp]
  if(length(rc)==0){
    sapply(1:(n+1),function(i) C_M[,i] <- x)
    R_L <- lapply(1:n,function (x) return(1:ncol(crn$mr)))
    return(list(C_M=C_M,R_L=R_L))
  }
  
  S <- crn$mp[sp,rc,drop=F]-crn$mr[sp,rc,drop=F]
  Sr <- crn$mr[sp,rc,drop=F]
  R_L <<- list()
  
  vp <- function(x_i){
    
    asp <- which(Sr%*%v > x_i)
    if(length(asp)>0){
      
    r_p <- x_i/Sr%*%v
    p <- rep(1,length(v))
    ind <- c()
    for (i in 1:ncol(S)){
      nsp <- which(Sr[,i]!=0)  
      if(length(nsp)!=0) p[i] <- max(0,min(r_p[nsp]))
      
    }
    R_L <- c(R_L,list(which(p==0)))
    vo <- v
    vo[p!=1] <- v[p!=1]*p[p!=1]
    return(vo)
  }
  else{return(v)}
  }
  for (i in 1:(n+1)){
    
    if(i==1){ C_M[,i] <- x; V_M[,i] <- v}
    else{ 
      C_M[,i] <- round(sapply(S%*%vp(C_M[,i-1]) + C_M[,i-1],function (x) max(0,x)),2)
      V_M[,i] <- round(vp(C_M[,i-1]),2)
      P_M[,i-1] <- round(S%*%vp(C_M[,i-1]),2)
    }
  }
  return(list(C_M=C_M,V_M=V_M,P_M=P_M))
}


rn.dina_step_stoi_2 <- function(crn,sp,x,v,ndp,n){
  
  # sp <- sp[x>0]
  rc <- c()
  for(i in 1:ncol(crn$mr)) if(all(crn$mp[-sp,i]==0) && all(crn$mr[-sp,i]==0)) rc <- c(rc,i) 
  C_M <- matrix(0,ncol=(n+1),nrow=length(x))
  V_M <- matrix(0,ncol=(n+1),nrow=length(v))
  R_L <- list()
  rownames(C_M) <- crn$sp.id[sp]
  if(length(rc)==0){
    sapply(1:(n+1),function(i) C_M[,i] <- x)
    R_L <- lapply(1:n,function (x) return(1:ncol(crn$mr)))
    return(list(C_M=C_M,R_L=R_L))
  }
  
  S <- crn$mp[sp,rc,drop=F]-crn$mr[sp,rc,drop=F]
  Sr <- crn$mr[sp,rc,drop=F]
  R_L <<- list()
  
  ndsp <- unique(unlist(sapply(ndp, function(x) which(crn$mr[,x]!=0))))
  
  vp <- function(x_i){
    
    x_nd <- (Sr%*%v)
    x_nd[setdiff(sp,ndsp)] <- 0
    r_p <- (x_i-x_nd)/(Sr%*%v-x_nd)
    r_p[(Sr%*%v-x_nd)==0] <- 1
    
    if(any(r_p[setdiff(sp,ndsp)]<1)){
      p <- rep(1,length(v))
      ind <- which(rc %in% ndp)
      for (i in setdiff(1:ncol(S),ind)){
        nsp <- which(Sr[,i]!=0)  
        p[i] <- max(0,min(1,r_p[nsp]))
        
      }
      
      R_L <- c(R_L,list(which(p==0)))
      vo <- v
      vo[p!=1] <- v[p!=1]*p[p!=1]
      return(vo)
    }
    else{return(v)}
  }
  for (i in 1:(n+1)){
    
    if(i==1){ C_M[,i] <- x; V_M[,i] <- v}
    else{ 
      C_M[,i] <- sapply(S%*%vp(C_M[,i-1]) + C_M[,i-1],function (x) max(0,x))
      V_M[,i] <- vp(C_M[,i-1])
    }
  }
  return(list(C_M=C_M,V_M=V_M,R_L=R_L))
  
}

# discrete dynamic concentration calculation function whit some of them governet by MAK
rn.dina_step_mix <- function(crn,sp,mak_r,x,v,n){
  
  # sp <- sp[x>0]
  rc <- c()
  
  for(i in 1:ncol(crn$mr)) if(all(crn$mp[-sp,i]==0) && all(crn$mr[-sp,i]==0)) rc <- c(rc,i) 
  C_M <- matrix(0,ncol=(n+1),nrow=length(x))
  V_M <- matrix(0,ncol=(n+1),nrow=length(v))
  R_L <- list()
  rownames(C_M) <- crn$sp.id[sp]
  if(length(rc)==0){
    sapply(1:(n+1),function(i) C_M[,i] <- x)
    R_L <- lapply(1:n,function (x) return(1:ncol(crn$mr)))
    return(list(C_M=C_M,R_L=R_L))
  }
  
  S <- crn$mp[sp,rc,drop=F]-crn$mr[sp,rc,drop=F]
  Sr <- crn$mr[sp,rc,drop=F]
  R_L <<- list()
  
  
  vp <- function(x){
    p <- v
    for(i in mak_r){
      p[i] <- 1
      if (length(crn$mr[,i])!=0) p[i] <- prod(x^crn$mr[,i])
      p[i] <- p[i]*v[i]
      # if(i==3) print(p[i])
    }
    # print("p:")
    # print(p[c(2,9,11)])
    return(p)
  }
  
  Sm<- function(x){
    asp <- which(x >= as.integer(round(Sr%*%vp(x))))
    # print(setdiff(sp,asp))
    rc <- c()
    #for(i in 1:ncol(crn$mr)) if(all(Sr[Sr[,i]!=0,i]*v[i]<x[Sr[,i]!=0] & x[Sr[,i]!=0]>0)) rc <- c(rc,i)
    for(i in 1:ncol(crn$mr)) if(all(S[-asp,i]==0) && all(crn$mr[-asp,i]==0)) rc <- c(rc,i)
    rc <- unique(union(rc,mak_r))
    St <- S
    R_L <<- c(R_L,list(setdiff(1:ncol(S),rc)))
    St[,setdiff(1:ncol(S),rc)] <-0
    return(St)
  }
  
  
  for (i in 1:(n+1)){
    if(i==1) C_M[,i] <- x
    else C_M[,i] <- sapply(Sm(C_M[,i-1])%*%vp(C_M[,i-1]) + C_M[,i-1], function(x) max(0,x))
    V_M[,i] <- vp(C_M[,i-1])
  }
  return(list(C_M=C_M,R_L=R_L,V_M=V_M))
  
}

# discrete dynamic concentration calculation function whit some of them governet by MAK
rn.dina_edo_mix <- function(crn,sp,mak_r,x,v,tf,dt){
  
  # sp <- sp[x>0]
  rc <- c()
  
  for(i in 1:ncol(crn$mr)) if(all(crn$mp[-sp,i]==0) && all(crn$mr[-sp,i]==0)) rc <- c(rc,i) 
  if(length(rc)==0){
    return(NULL)
  }
  
  S <- crn$mp[sp,rc,drop=F]-crn$mr[sp,rc,drop=F]
  Sr <- crn$mr[sp,rc,drop=F]
  
  
  vp <- function(x){
    p <- v
    for(i in mak_r){
      p[i] <- 1
      if (length(crn$mr[,i])!=0) p[i] <- prod(x^crn$mr[,i])
      p[i] <- p[i]*v[i]
      # if(i==3) print(p[i])
    }
    # print("p:")
    # print(p[2,9,11])
    # return(p)
  }
  
  Sm<- function(x){
    asp <- which(x >= as.integer(round(Sr%*%vp(x))))
    # print(setdiff(sp,asp))
    rc <- c()
    #for(i in 1:ncol(crn$mr)) if(all(Sr[Sr[,i]!=0,i]*v[i]<x[Sr[,i]!=0] & x[Sr[,i]!=0]>0)) rc <- c(rc,i)
    for(i in 1:ncol(crn$mr)) if(all(S[-asp,i]==0) && all(crn$mr[-asp,i]==0)) rc <- c(rc,i)
    rc <- unique(union(rc,mak_r))
    St <- S
    # R_L <<- c(R_L,list(setdiff(1:ncol(S),rc)))
    St[,setdiff(1:ncol(S),rc)] <-0
    return(St)
  }
  
  mix_dina <- function(t, y, parms){
    dy <- Sm(y)%*%vp(y)
    return(list(dy))
  }
  
  times <- seq(from = 0, to = tf, by = dt)
  C_M  <- ode(times = times, y = x, func = mix_dina, parms = NULL)
  
  return(C_M)
}

rn.is_ssm <- function(crn=rn, sp=crn$sp.id){
  
  reac_r <- which(sapply(1:length(crn$reac),function(r) all(crn$reac[[r]]$reactants %in% sp)))
  reac_sp <- unique(unlist(sapply(reac_r, function(x) crn$reac[[x]]$reactants)))
  prod_sp <- unique(unlist(sapply(reac_r, function(x) crn$reac[[x]]$products)))
  
  if(length(setdiff(reac_sp,prod_sp))==0) return(T)
  else return(F)
}