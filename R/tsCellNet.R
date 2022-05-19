library(igraph)
library(homologene)
library(lattice)
library(Hmisc)
library(ggpubr)
library(pheatmap)
library(reshape2)
library(ggalluvial)
library(tseries)
#' tsCellNet
#'
#' construction time-series cell-cell interactions using scRNA-seq data
#'
#' @param datloc three input file path
#' @param expnam expression file name
#' @param timnam time-point file name
#' @param typnam cell-type file name
#' @param resloc result store path
#' @param mainlab result prefix name
#' @param species species name
#' @param controltime the starting point name
#' @param treattime if or not need to compare cell-cell interactions with starting point name. If not set the value as "total".

#' @export
#' @importFrom grDevices colorRampPalette dev.off pdf
#' @importFrom stats aggregate na.omit sd time
#' @importFrom utils head read.csv write.csv
#' @examples
#' out<-datainput(datloc,expnam,typnam,timnam,"FALSE",resloc)
#' inpexp<-as.data.frame(out$exp) #[,kep]
#' inptyp<-as.data.frame(out$typ)
#' colnames(inptyp)<-colnames(out$typ)
#' rownames(inptyp)<-rownames(out$typ) 
#' inptim<-as.data.frame(out$tim) 
#' usedcolo<-definecolo(tsccnetcol,inptyp$Type,resloc,mainlab)
#' sankey(inptim$Time,inptyp$Type,usedcolo,mainlab,resloc)
#' memsercpros<-preparepairs(inpexp,species,dbloc,mainlab)
#' singlenet(inpexp,inptyp,"Type",inptim,"Time",memsercpros,usedcolo,"FALSE",dbloc,mainlab) 
#' developcorrelation(inptim,usedcolo,resloc)
#' comparednet(controltime,treattime,usedcolo,resloc)
#' importance_plot(resloc,"row","outgoing","None")

dbloc="./database/"
load(paste(dbloc,"tsCCNETcolor.RData",sep=""))
load(paste(dbloc,"membrane_protein_mouse.RData",sep=""))
load(paste(dbloc,"secreted_protein_mouse.RData",sep=""))

normalized <- function(x){
  y<-(x-min(x))/(max(x)-min(x))
  return(y)
}

convertGeneList <- function(x,sourcespecies,aimspecies){
  speciesid<-homologene::taxData$tax_id[homologene::taxData$name_txt==sourcespecies]
  aimspeciesid<-homologene::taxData$tax_id[homologene::taxData$name_txt==aimspecies]
  memout<-homologene(x,inTax = speciesid,outTax = aimspeciesid)
  genesV2<-na.omit(memout[,1:2])
  
  humanx <- genesV2[, 2]
  names(humanx)<-genesV2[,1]
  print(head(humanx))
  return(humanx)
}

datainput <- function(datloc,expnam,typnam,timnam,ifQC,resloc){
  shell_cmd<-paste("mkdir ",resloc,sep = "")
  grep_out<-system(shell_cmd, intern = TRUE)
  cat(grep_out)
  shell_cmd<-paste("mkdir ",resloc,"tmp",sep = "")
  grep_out<-system(shell_cmd, intern = TRUE)
  cat(grep_out)
  tmploc=paste(resloc,"tmp",sep = "")
  datexp<-read.csv(paste(datloc,expnam,sep=""),check.names = F,header = T,row.names = 1)
  dattyp<-read.csv(paste(datloc,typnam,sep=""),check.names = F,header = T,row.names = 1)
  dattim<-read.csv(paste(datloc,timnam,sep=""),check.names = F,header = T,row.names = 1)
  out<-list()
  if(ifQC=="True"|ifQC=="TRUE"){
  colnamytyp<-colnames(ytyp)
  colnamytim<-colnames(ytim)
  sumCols<-colSums(y)
  print(paste("Quality control:dimension of data(raw) is ",dim(y)[1],dim(y)[2],sep=" "))
  y<-y[,sumCols>5000]
  ytyp<-as.data.frame(ytyp[sumCols>5000,])
  ytim<-as.data.frame(ytim[sumCols>5000,])
  print(paste("Quality control:dimension of data(at least 5K reads) is ",dim(y)[1],dim(y)[2],sep=" "))
  seqdepth<-colSums(y)/10^6
  y<-sweep(y,2,seqdepth,"/")
  print("Quality control:remove sequence depth among cells")
  genenum<-apply(y,1,function(x) sum(x>0))
  y<-y[genenum>5,]
  print(paste("Quality control:dimension of data(at least express 5 cells) is ",dim(y)[1],dim(y)[2],sep=" "))
  celnum<-apply(y,2,function(x) sum(x>0))
  y<-y[,celnum>500]
  ytyp<-as.data.frame(ytyp[celnum>500,])
  ytim<-as.data.frame(ytim[celnum>500,])
  colnames(ytyp)<-colnamytyp
  rownames(ytyp)<-colnames(y)
  colnames(ytim)<-colnamytim
  rownames(ytim)<-colnames(y)
  print(paste("Quality control:dimension of data(at least express 500 genes) is ",dim(y)[1],dim(y)[2],sep=" "))
  out$exp<-y
  out$typ<-ytyp
  out$tim<-ytim
  return(out)
  }else if(ifQC=="False"|ifQC=="FALSE"){
  out$exp<-datexp
  out$typ<-dattyp
  out$tim<-dattim
  return(out)
  }
  print("Finish quality control!") 
}

memserchomotrans<-function(mempros,sercpros,aimspecies){
  sourcespecies<-"Mus musculus"
  aimnam<-substr(aimspecies, 1, 3)
  memsercpros<-list()
  if(aimspecies!="Mus musculus" & aimspecies %in% homologene::taxData$name_txt){
    memout<-convertGeneList(mempros,sourcespecies,aimspecies)
    memprot<-na.omit(memout)
    sercout<-convertGeneList(sercpros,sourcespecies,aimspecies)
    sercprot<-na.omit(sercout)
    memsercpros$memprot<-memprot
    memsercpros$sercprot<-sercprot
    print(paste("Step0:Number of used membrane/secreted proteins of",aimnam,"is ",length(memprot),length(sercprot),sep=" "))
    return(memsercpros)
  }else if(aimspecies=="Mus musculus"){
    memprot<-mempros
    sercprot<-sercpros 
    memsercpros$memprot<-memprot
    memsercpros$sercprot<-sercprot
    print(paste("Step0:Number of used membrane/secreted proteins of",aimnam,"is ",length(memprot),length(sercprot),sep=" "))
    return(memsercpros)
  }else{
    print("Error1:check the species dataset name or don't support this species! Avaliable species:") 
    print(homologene::taxData$name_txt)
    break  
  }
}

LRdbhomotrans<-function(yexp,mempros,sercpros,aimspecies,dbloc){
  sourcespecies="Homo sapiens"
  aimnam<-substr(aimspecies, 1, 3)
  if(aimspecies!="Homo sapiens" & aimspecies %in% homologene::taxData$name_txt){
    memprot<-mempros
    sercprot<-sercpros
    load(paste(dbloc,"LRdb/lr-tsCC-human.Rdata",sep=""))
    pairsdata<-tsCCPairsDb
    genel<-as.character(unlist(pairsdata[,1]))
    gener<-as.character(unlist(pairsdata[,2]))
    genelout<-convertGeneList(genel,sourcespecies,aimspecies)
    generout<-convertGeneList(gener,sourcespecies,aimspecies)
    
    res<-tsCCPairsDb
    res$Ligand<-genelout[pairsdata[,1]]
    res$Receptor<-generout[pairsdata[,2]]
    usedpairsdata0<-na.omit(res)
    usedpairsdata1<-usedpairsdata0[(usedpairsdata0[,1] %in% rownames(yexp))  & (usedpairsdata0[,2] %in% rownames(yexp)),]
    usedpairsdata<-usedpairsdata1[(usedpairsdata1[,1] %in% unique(c(memprot,sercprot)))  & (usedpairsdata1[,2] %in% unique(c(memprot,sercprot))),]
    return(usedpairsdata)
    print(paste("Step1:dimension of used L&R pairs of",aimnam,"is ",dim(usedpairsdata)[1],sep=" "))
    print(paste("Step1:dimension of used L R genes of",aimnam,"is ",length(unique(usedpairsdata[,1])),length(unique(usedpairsdata[,2])),sep=" "))
    
  }else if(aimspecies %in% c("Homo sapiens","Mus musculus")){
    if(aimspecies=="Homo sapiens"){
    load(paste(dbloc,"LRdb/lr-tsCC-human.Rdata",sep=""))
    }else if(aimspecies=="Mus musculus"){
    load(paste(dbloc,"LRdb/lr-tsCC-mouse.Rdata",sep=""))  
    }
    pairsdata<-tsCCPairsDb
    
    memprot<-mempros
    sercprot<-sercpros 
    usedpairsdata0<-pairsdata
    print(paste("Step1:dimension of used L&R pairs of",aimnam,"is ",dim(usedpairsdata0)[1],sep=" "))
    usedpairsdata1<-usedpairsdata0[(usedpairsdata0[,1] %in% rownames(yexp))  & (usedpairsdata0[,2] %in% rownames(yexp)),]
    usedpairsdata<-usedpairsdata1[(usedpairsdata1[,1] %in% unique(c(memprot,sercprot)))  & (usedpairsdata1[,2] %in% unique(memprot)),]
    return(usedpairsdata)
    print(paste("Step1:dimension of expressed L&R pairs of",aimnam,"is ",dim(usedpairsdata)[1],sep=" "))
    print(paste("Step1:dimension of expressed L R genes of",aimnam,"is ",length(unique(usedpairsdata[,1])),length(unique(usedpairsdata[,2])),sep=" "))
    
  }else{
    print("Error1:missing this species! Avaliable species:") 
    print(homologene::taxData$name_txt)
    break  
  }
}

select_base_exp<-function(inputs){
  n=0
  inputs$Group.1<-as.character(inputs$Group.1) 
  inputs$Group.2<-as.character(inputs$Group.2) 
  inputs$variable<-as.character(inputs$variable) 
  inputs$value<-as.numeric(inputs$value)
  colnames(inputs)<-c("Group.1","Group.2","variable","value")
  outputs<-as.data.frame(matrix(NA,nc=4,nr=1))
  line1<-0
  L1<-inputs
  LL1<-L1
  L1line<-mean(LL1$value)
  s0<-0
  got=0
  if(max(LL1$value)>=L1line & got==0 & L1line!=0){ 
    for(n in 1:length(which(LL1$value==max(LL1$value)))){
    line1<-line1+1
    outputs[line1,]=LL1[which(LL1$value==max(LL1$value))[n],] 
    got=1
    }
  }else if(max(LL1$value)>=L1line & got==0 & L1line==0){ 

  }else{ 
    for(i in 1:dim(LL1)[1]){
        if(s0<L1line & got==0){ 
          line1<-line1+1
          outputs[line1,]=LL1[which(LL1$value==max(LL1$value)),] 
          L1l=which(LL1$value==max(LL1$value))
          s0<-s0+max(LL1$value) 
          LL1<-LL1[-L1l,] 
        }else if(s0>=L1line & got==0){
          line1<-line1+1
          outputs[line1,]=LL1[which(LL1$value==max(LL1$value)),] 
          got=got+1
          }
        }
  }
  return(outputs)
}

select_base_pval<-function(totalLRone,inputL,inputR,lenL,lenR){ 
set.seed(1)
pval<-c()
permres<-c()
for(dd in 1:100){ 
    idl=sample(1:ncol(inputL),lenL,replace =T)
    idr=sample(1:ncol(inputR),lenR,replace =T)
    permres[dd]<-(sum(unlist(inputL[idl]))*sum(unlist(inputR[idr]))) 
    }
    pval=(sum(permres>(totalLRone$Lexp*totalLRone$Rexp))/length(permres))
    return(pval)
}

single_network_plot<-function(kepLRpvalinp,mainlab,usedcolo){
library(igraph)
pdf(paste(resloc,mainlab,"net.pdf",sep=""))
  
  set.seed(200)
  plotinp0<-kepLRpvalinp[,colnames(kepLRpvalinp) %in% c("Lloc","Rloc")]
  plotinp0$val<-1
  plotinp<-aggregate(plotinp0$val,by=list(plotinp0$Lloc,plotinp0$Rloc),sum) 
  colnames(plotinp)<-c("from","to","edge")
  write.csv(plotinp,paste(resloc,mainlab,"netinput.csv",sep=""),quote = F)
  
  g <- graph_from_data_frame(plotinp, directed=TRUE)
  
  usedcolo1<-usedcolo[unique(c(plotinp$from,plotinp$to))]
  
  V(g)$color <-usedcolo1
  E(g)$width <-normalized(E(g)$edge)*10+1
  E(g)$arrow.size <- 0.3
  E(g)$edge.color <- "gray80"
  l <-do.call("layout_in_circle", list(g)) 
  plot(g,layout=l,margin=0.3,vertex.label.dist=2,edge.curved=0.5,vertex.frame.color=NA,main=mainlab)
  inlink<-as.data.frame(strength(g, mode="in",weights = E(g)$edge))
  outlink<-as.data.frame(strength(g, mode="out",weights = E(g)$edge))
  totallink<-as.data.frame(strength(g, mode="total",weights = E(g)$edge))
  betweennesscore<-as.data.frame(betweenness(g, directed=T, weights=E(g)$edge))
  strengout<-as.data.frame(cbind(inlink,outlink,totallink,betweennesscore))
  colnames(strengout)<-c("incoming","outgoing","total","betweennesscore")
  rownames(strengout)<-paste(mainlab,rownames(strengout),sep = ":")
  write.csv(strengout,paste(resloc,mainlab,"netdetail.csv",sep=""),quote = F)
  
  hs <- hub_score(g, weights=E(g)$edge)$vector
  plot(g,margin=0.3,vertex.label.dist=2,edge.curved=0.5,layout=l,vertex.frame.color=NA, vertex.size=hs*20, main=paste(mainlab,"Hub net"))
dev.off()
}

single_network_gene_plot<-function(kepLRpvalinp,usedgene,mainlab,usedcolo){
  library(igraph)
  pdf(paste(resloc,usedgene,mainlab,"-net.pdf",sep=""))
  
  set.seed(200)
  kepLRpvalinp1<-kepLRpvalinp
  inp0<-kepLRpvalinp1[(kepLRpvalinp$Lgen %in% usedgene | kepLRpvalinp$Rgen %in% usedgene),]
  plotinp0<-as.data.frame(kepLRpvalinp1[(kepLRpvalinp$Lgen %in% usedgene | kepLRpvalinp$Rgen %in% usedgene),colnames(kepLRpvalinp) %in% c("Lloc","Rloc")])
  plotinp0$LRgene<-paste(inp0$Lgen,inp0$Rgen,sep=":")
  plotinp0$val<-1
  plotinp<-aggregate(plotinp0$val,by=list(plotinp0$Lloc,plotinp0$Rloc),sum) 
  for(s in 1:dim(plotinp)[1]){
    plotinp[s,4]=paste(plotinp0$LRgene[plotinp0$Lloc==plotinp$from[s] & plotinp0$Rloc==plotinp$to[s]],collapse = "/")
  }
  colnames(plotinp)<-c("from","to","edge","label")
  g <- graph_from_data_frame(plotinp[,1:3], directed=TRUE)
  
  usedcolo1<-usedcolo[unique(c(plotinp$from,plotinp$to))]
  
  V(g)$color <-usedcolo1
  E(g)$width <-normalized(E(g)$edge)*10+1
  E(g)$arrow.size <- 0.3
  E(g)$edge.color <- "gray80"
  l <-do.call("layout_in_circle", list(g))
  plot(g,layout=l,margin=0.3,vertex.label.dist=2,edge.curved=0.5,vertex.frame.color=NA,main=mainlab)
  write.csv(plotinp,paste(resloc,usedgene,mainlab,"-netdetail.csv",sep=""),quote = F)
  
  hs <- hub_score(g, weights=E(g)$edge)$vector
  plot(g,margin=0.3,vertex.label.dist=2,edge.curved=0.5,layout=l,vertex.frame.color=NA, vertex.size=hs*20, main=paste(mainlab,"Hub net"))
  dev.off()
}

stat_LRcor<-function(dataL,dataR,dattim){
  library(Hmisc)
  d=0
  tim<-as.character(unique(dattim$Time[order(dattim$Order)]))
  res<-as.data.frame(matrix(NA,nr=1,nc=5))
  colnames(res)<-c("Ligand","Receptor","p","r","L:R")
  for(i in 1:length(unique(dataL$variable))){
    inp1n<-as.data.frame(dataL[as.character(dataL$variable)==unique(dataL$variable)[i] ,])
    inp2n<-as.data.frame(dataR[as.character(dataL$variable)==unique(dataL$variable)[i] ,])
    
    inp1<-reshape2::dcast(inp1n,Group.2~Group.1)
    inp1s<-inp1[,-1]
    rownames(inp1s)<-inp1[,1]
    inp1s[is.na(inp1s)]=0
    inp1sord<-inp1s[tim,]
    inp2<-reshape2::dcast(inp2n,Group.2~Group.1) 
    inp2s<-inp2[,-1]
    rownames(inp2s)<-inp2[,1]
    inp2s[is.na(inp2s)]=0
    inp2sord<-inp2s[tim,]
    
    insd1 = apply(inp1sord, 2, sd, na.rm = TRUE)
    insd2 = apply(inp2sord, 2, sd, na.rm = TRUE)
    
    inp1sd<-as.data.frame(inp1sord[,insd1>0.2])
    inp2sd<-as.data.frame(inp2sord[,insd2>0.2])
    
    if(dim(inp1sd)[2]==0 | dim(inp2sd)[2]==0){
    }else{
    
    rownames(inp1sd)<-rownames(inp1sord)
    colnames(inp1sd)<-colnames(inp1sord)[insd1>0.2]
    rownames(inp2sd)<-rownames(inp2sord)
    colnames(inp2sd)<-colnames(inp1sord)[insd2>0.2]
    
    if(dim(inp1sd)[2]>2 & dim(inp2sd)[2]>2){ 
    cores0<-rcorr(as.matrix(inp1sd),as.matrix(inp2sd))
    coresr<-as.data.frame(cores0$r)[1:ncol(inp1sd),(ncol(inp1sd)+1):ncol(cores0$r)]
    coresP<-as.data.frame(cores0$P[1:ncol(inp1sd),(ncol(inp1sd)+1):ncol(cores0$r)])
    comnam<-paste(rep(colnames(inp1sd),each=length(colnames(inp2sd))),rep(colnames(inp2sd),length(colnames(inp1sd))),sep=":")
    coresr[is.na(coresr)]=0
    coresP[is.na(coresP)]=1
    if(length(coresP[coresP<0.01])>0){
    for(w in 1:length(coresP[coresP<0.01])){
      d=d+1
      res[d,1]=sub("\\.[0-9]+$","",as.character(unique(dataL$variable)[i])) 
      res[d,2]=sub("\\.[0-9]+$","",as.character(unique(dataR$variable)[i])) 
      res[d,3]=round(coresP[coresP<0.01][w],4) 
      res[d,4]=round(coresr[coresP<0.01][w],4) 
      res[d,5]=comnam[coresP<0.01][w] 
    }
    }
    }
    }
  }
  return(na.omit(res))
  write.csv(na.omit(res),paste(resloc,"significant_LRcor.csv",sep=""),quote = F,row.names = F)
}

plot_LR_exps_label<-function(dataL,dataR,sigres,dattim,usedcolo){ 
  library(ggpubr)
  pdf(paste(resloc,"plot_LR_exps_line.pdf",sep = ""))
  for(i in 1:length(unique(dataL$variable))){ 
    inp1n<-as.data.frame(dataL[as.character(dataL$variable)==as.character(unique(dataL$variable)[i]) ,])
    tim<-as.character(unique(dattim$Time[order(dattim$Order)]))
    timord<-as.numeric(unique(dattim$Order[order(dattim$Order)]))
    inp1n$Group.2<-factor(inp1n$Group.2,levels =tim)
    inp1n$Group.2<-as.numeric(inp1n$Group.2)
    siginp<-sigres[sigres$Ligand==as.character(unique(dataL$variable)[i]) & sigres$Receptor==as.character(unique(dataR$variable)[i]),]
    inp1n$variable<-as.character(inp1n$variable)
    usedcolo1<-usedcolo[unique(inp1n$Group.1)] 
    p1<-ggplot(data=inp1n,mapping=aes(x=Group.2,y=value,colour = Group.1))+ geom_point()+theme_classic()+ geom_line()+ scale_x_continuous(breaks=timord, labels=tim)+labs(title=sub("\\.[0-9]+$","",as.character(unique(dataL$variable)[i])),y="Expression of Ligand",x="Time")+theme(axis.title.y =element_text(size=10),axis.text.x  =element_text(angle = 15),legend.position = "none")+scale_color_manual(values=usedcolo1,name="Type") 
    
    if(dim(siginp)[1]==0){
    p2<-NULL
    }else{
    p2<-ggtexttable(siginp, rows = NULL, theme = ttheme("light"))
    }
    
    inp2n<-as.data.frame(dataR[as.character(dataR$variable)==as.character(unique(dataR$variable)[i]) ,])
    inp2n$Group.2<-factor(inp2n$Group.2,levels =tim)
    inp2n$Group.2<-as.numeric(inp2n$Group.2)
    usedcolo1<-usedcolo[unique(inp2n$Group.1)]
    p3<-ggplot(data=inp2n,mapping=aes(x=Group.2,y=value,colour = Group.1))+ geom_point()+theme_classic()+ geom_line()+ scale_x_continuous(breaks=timord, labels=tim)+labs(title=sub("\\.[0-9]+$","",as.character(unique(dataR$variable)[i])),y="Expression of Receptor",x="Time")+theme(axis.title.y =element_text(size=10),legend.position = "none")+scale_color_manual(values=usedcolo1,name="Type")+theme(axis.text.x  =element_text(angle = 15))
    
    p<-ggarrange(p2,
              ggarrange(p1, p3, ncol = 2, labels = c("B", "C")),
              nrow = 2,
              labels = "A", heights  = c(0.8,1.9)) 
    
    print(p)
    if(i==1|i==2|i==3){
    usedcolo1<-usedcolo[unique(inp2n$Group.1)]  
    p0<-ggplot(data=inp2n,mapping=aes(x=Group.2,y=log(value),colour = Group.1))+ geom_point()+theme_classic()+ geom_line()+ scale_x_continuous(breaks=timord, labels=tim)+labs(title=sub("\\.[0-9]+$","",as.character(unique(dataR$variable)[i])),y="Log transformed expression of Receptor",x="Time")+theme(axis.title.y =element_text(size=10))+scale_color_manual(values=usedcolo1,name="Type")
    print(p0)  
    }
  }
  dev.off()
}

compare_network_plot<-function(kepLRpval,controltime,usedcolo){
  controlinp<-kepLRpval[kepLRpval$time==controltime,]
  otherinp<-kepLRpval[kepLRpval$time!=controltime,]
  for(t in 1:length(unique(otherinp$time))){
    tmpinp<-otherinp[otherinp$time==unique(otherinp$time)[t],]
    incinp<-tmpinp[!((tmpinp$Lloc %in% controlinp$Lloc) & (tmpinp$Rloc %in% controlinp$Rloc) & (tmpinp$Lgen %in% controlinp$Lgen) & (tmpinp$Rgen %in% controlinp$Rgen)),]
    write.csv(incinp,paste(resloc,unique(otherinp$time)[t],"vs",controltime,"-comparenet-increased.csv",sep=""),quote = F,row.names = F)
    decinp<-controlinp[!((controlinp$Lloc %in% tmpinp$Lloc) & (controlinp$Rloc %in% tmpinp$Rloc) & (controlinp$Lgen %in% tmpinp$Lgen) & (controlinp$Rgen %in% tmpinp$Rgen)),]
    write.csv(decinp,paste(resloc,unique(otherinp$time)[t],"vs",controltime,"-comparenet-decreased.csv",sep=""),quote = F,row.names = F)
    
    single_network_plot(incinp,paste("Compared with",controltime,"network in",unique(otherinp$time)[t]),usedcolo)
  }
}

dot_exp_pval<-function(kepLRpval,dattimpval){
pdf(paste(resloc,"plot_LR_exps_dot.pdf",sep = ""))
  dotinp<-kepLRpval[,colnames(kepLRpval) %in% c("time","pval")]
  dotinp$LRgene<-paste(kepLRpval$Lgen,kepLRpval$Rgen,sep=":")
  dotinp$LRloc<-paste(kepLRpval$Lloc,kepLRpval$Rloc,sep=":")
  dotinp$LRexp<-kepLRpval$Lexp*kepLRpval$Rexp 
  tim<-as.character(unique(dattimpval$Time[order(dattimpval$Order)]))
  timord<-unique(dattimpval$Order[order(dattimpval$Order)])
  dotinp$time<-factor(dotinp$time,levels =tim)
  
  p1<-ggplot(dotinp, aes(time, LRgene)) + geom_point(aes(size=-log10(dotinp$LRexp+0.000001),color=dotinp$pval))+scale_size_continuous(range=c(0.1,3))+ facet_wrap( ~ dotinp$LRloc)+ theme(axis.text.y = element_text(size = 7),strip.text = element_text(size = 7),axis.text.x=element_text(face="bold",size=8,angle=-45))+theme(strip.background = element_blank(), strip.placement = "outside",legend.title=element_text("pvalue"))+labs(x="Time",y="Ligand:Receptor",color="pval",size="Exp")
  p1
dev.off()
}

random_type_treat<-function(datexp,dattim,dattyp){
  for(t in 1:length(unique(dattim$Time))){
    for(c in 1:length(unique(dattyp$Type))){
      in0<-datexp[,dattim$Time==unique(dattim$Time)[t] & dattyp$Type==unique(dattyp$Type)[c]]
      for(i in 1:10){
      a <- rowSums(inn0[, sample(1:ncol(inn0), round(ncol(inn0)*0.5))], na.rm = T)
      if(i == 1){
        id_tmp <- data.frame(exp = a)
      }else{
        id_tmp <- cbind.data.frame(id_tmp, as.data.frame(a))
      }
      }
      colnames(id_tmp) <- paste0(as.character(unique(dattim$Time)[t]),":",as.character(unique(dattyp$Type)[c]),":", 1:100)
    }
      if(t==1){
        totdat<-id_tmp
      }else{
        totdat<-cbind.data.frame(totdat,id_tmp)
      }
  }
  return(totdat)
  }

heat_cellLR<-function(totdat,celltypeL,celltypeR,pairs,dattim){
  types<-sapply(strsplit(colnames(totdat)),"[[",2)
  times<-sapply(strsplit(colnames(totdat)),"[[",1)
  heatinpL0<-totdat[rownames(totdat)%in% unique(unlist(pairs[,1])), types==celltypeL]
  dattimL0<-times[types==celltypeL]
  heatinpR0<-totdat[rownames(totdat)%in% unique(unlist(pairs[,2])), types==celltypeR]
  dattimR0<-times[types==celltypeR]
  
  weekL<-dattimL0
  typeL<-celltypeL
  weekcol<-brewer.pal(9,"YlGn")
  names(weekcol)<-as.character(unique(dattim$Time[order(dattim$Order)]))
  typecol<-tsccnetcol
  names(typecol)<-as.character(unique(celltypeL))
  
  annotation_colL = data.frame(TypeClass = typeL,Week=weekL)
  rownames(annotation_col) = colnames(heatinpL0)
  ann_colors = list(TypeClass = typecol,Week = weekcol)
  
  pl<-pheatmap(heatinpL0,annotation_colors = ann_colors,annotation_col=annotation_colL,show_rownames=T,show_colnames=F,scale = "row",cluster_cols = F)
  
  weekR<-dattimR0
  typeR<-celltypeR
  annotation_colR = data.frame(TypeClass = typeR,Week=weekR)
  rownames(annotation_col) = colnames(heatinpR0)
  pr<-pheatmap(heatinpR0,annotation_colors = ann_colors,annotation_col=annotation_colR,show_rownames=T,show_colnames=F,scale = "row",cluster_cols = F)
  library(gridExtra)
  grid.arrange(pl, pr, nrow = 1)
}

preparepairs<-function(yexp,aimspecies,dbloc,mainlab){
  memsercpros<-memserchomotrans(mempros,sercpros,aimspecies)
  print(paste("Membrane protein num is ",length(unique(memsercpros$memprot)),sep=""))
  print(paste("Secreted protein num is ",length(unique(memsercpros$sercprot)),sep=""))
  usedpairsdata<-LRdbhomotrans(yexp,memsercpros$memprot,memsercpros$sercprot,aimspecies,dbloc)
  print(paste("Total L&R pairs is ",dim(usedpairsdata)[1],"; Total Ligand genes ",length(unique(usedpairsdata[,1])),"; Total Receptor genes ",length(unique(usedpairsdata[,2])),sep=" "))
  write.csv(usedpairsdata,paste(dbloc,mainlab,"-usedpairs.csv",sep=""),quote = F)
  return(memsercpros)
  }

singlenet<-function(yexp,ytyp,ytypnam,ytim,ytimnam,memsercpros,usedcolo,scaled,dbloc,mainlab){
  usedpairsdata<-read.csv(paste(dbloc,mainlab,"-usedpairs.csv",sep=""),header = T,row.names = 1)
  print(paste("Total L&R pairs is ",dim(usedpairsdata)[1],"; Total Ligand genes ",length(unique(usedpairsdata[,1])),"; Total Receptor genes ",length(unique(usedpairsdata[,2])),sep=" "))
  
  lexp<-na.omit(yexp[usedpairsdata[,1],]) 
  rexp<-na.omit(yexp[usedpairsdata[,2],])
  print(paste("a.Used L&R pairs is ",dim(lexp)[1],sep=" "))
  if(scaled=="TRUE" | scaled=="True"){
    keepl<-rowSums(lexp > 0) >=10  
    keepr<-rowSums(rexp >0) >=10 
    lexps<-log((lexp[keepl & keepr,]+1),base=2) 
    rexps<-log((rexp[keepr & keepl,]+1),base=2)
  }else if(scaled=="FALSE" | scaled=="False"){
    lexps<-lexp
    rexps<-rexp  
  }
  print(paste("Used L&R pairs is ",dim(lexps)[1],"; Total cell num is",dim(lexps)[2],sep=" "))
  
  lexpsT0<-as.data.frame(t(lexps))
  lexpsT0$id=rownames(lexpsT0)
  rexpsT0<-as.data.frame(t(rexps))
  rexpsT0$id=rownames(rexpsT0)
  ytyp$id=rownames(ytyp)
  ytim$id=rownames(ytim)
  
  lexpsT1<-merge(lexpsT0,ytyp,by="id")
  lexpsT2<-merge(lexpsT1,ytim,by="id") 
  lexpsT<-lexpsT2[,-1] 
  rownames(lexpsT)<-lexpsT2[,1]
  rexpsT1<-merge(rexpsT0,ytyp,by="id")
  rexpsT2<-merge(rexpsT1,ytim,by="id")
  rexpsT<-rexpsT2[,-1]
  rownames(rexpsT)<-rexpsT2[,1]
  lg<-sub("\\.[0-9]+","",colnames(lexpsT))
  lgorder<-colnames(lexpsT)[!colnames(lexpsT) %in% c(ytypnam,ytimnam,"Order") ]
  inl1<-lexpsT[,lg %in% c(unique(memsercpros$memprot),ytypnam,ytimnam,"Order")] 
  inl2<-lexpsT[,!lg %in% c(unique(memsercpros$memprot))] 
  lexpsTa1<-aggregate(inl1[,!colnames(inl1) %in% c(ytypnam,ytimnam,"Order")],by=list(inl1[,colnames(inl1)==ytypnam],inl1[,colnames(inl1)==ytimnam]),FUN=mean)
  lexpsTa2<-aggregate(inl2[,!colnames(inl2) %in% c(ytypnam,ytimnam,"Order")],by=list(inl2[,colnames(inl2)==ytypnam],inl2[,colnames(inl2)==ytimnam]),FUN=sum)
  all(lexpsTa1$Group.1==lexpsTa2$Group.1)
  all(lexpsTa1$Group.2==lexpsTa2$Group.2)
  lexpsTam0<-cbind(lexpsTa1,lexpsTa2[,-c(1,2)])
  lexpsTam1<-lexpsTam0[,c("Group.1","Group.2",lgorder)]
  lexpsTa<-lexpsTam1
  data1<-melt(lexpsTa,id.var=c("Group.1","Group.2")) 
  rg<-sub("\\.[0-9]+","",colnames(rexpsT))
  rgorder<-colnames(rexpsT)[!colnames(rexpsT) %in% c(ytypnam,ytimnam,"Order") ]
  inr1<-rexpsT[,rg %in% c(unique(memsercpros$memprot),ytypnam,ytimnam,"Order")] 
  inr2<-rexpsT[,!rg %in% c(unique(memsercpros$memprot))] 
  rexpsTa1<-aggregate(inr1[,!colnames(inr1) %in% c(ytypnam,ytimnam,"Order")],by=list(inr1[,colnames(inr1)==ytypnam],inr1[,colnames(inr1)==ytimnam]),FUN=mean) 
  rexpsTa2<-aggregate(inr2[,!colnames(inr2) %in% c(ytypnam,ytimnam,"Order")],by=list(inr2[,colnames(inr2)==ytypnam],inr2[,colnames(inr2)==ytimnam]),FUN=sum) 
  all(rexpsTa1$Group.1==rexpsTa2$Group.1)
  all(rexpsTa1$Group.2==rexpsTa2$Group.2)
  rexpsTam0<-cbind(rexpsTa1,rexpsTa2[,-c(1,2)])
  rexpsTam1<-rexpsTam0[,c("Group.1","Group.2",rgorder)]
  rexpsTa<-rexpsTam1
  data2<-melt(rexpsTa,id.var=c("Group.1","Group.2"))
  write.csv(data1,paste(resloc,"tmp/","data1.csv",sep=""),quote = F,row.names = F)
  write.csv(data2,paste(resloc,"tmp/","data2.csv",sep=""),quote = F,row.names = F)
  n=0
  totalLR<-as.data.frame(matrix(NA,nc=8,nr=1))
  colnames(totalLR)<-c("time","Lloc","Lgen","Lexp","Rloc","Rgen","Rexp","type")
  for(t in 1:length(unique(data1$Group.2))){ 
    for(g in 1:length(unique(data1$variable))){ 
      Ltmp<-na.omit(select_base_exp(data1[(data1$variable==as.character(unique(data1$variable)[g])) & (data1$Group.2==as.character(unique(data1$Group.2)[t])),]))
      Rtmp<-na.omit(select_base_exp(data2[(data2$variable==as.character(unique(data2$variable)[g])) & (data2$Group.2==as.character(unique(data2$Group.2)[t])),]))
      colnames(Ltmp)<-colnames(data1)
      colnames(Rtmp)<-colnames(data2)
      if(dim(Rtmp)[1]!=0 & dim(Ltmp)[1]!=0){
        for(ll in 1: dim(Ltmp)[1]){
          for(rr in 1: dim(Rtmp)[1]){
            n=n+1
            totalLR[n,1]=as.character(unique(data1$Group.2)[t]) 
            totalLR[n,2]=Ltmp$Group.1[ll] 
            lg<-sapply(strsplit(Ltmp$variable[ll],".",fixed = T),"[[",1) 
            totalLR[n,3]=lg 
            totalLR[n,4]=Ltmp$value[ll] 
            totalLR[n,5]=Rtmp$Group.1[rr] 
            rg<-sapply(strsplit(Rtmp$variable[rr],".",fixed = T),"[[",1)
            totalLR[n,6]=rg 
            totalLR[n,7]=Rtmp$value[rr] 
            totalLR[n,8]=paste(Ltmp$Group.1[ll],Rtmp$Group.1[rr],sep=":") 
          }
        }
      }
    }
  }
  n=0
  totalLRpval<-as.data.frame(matrix(NA,nc=9,nr=1))
  colnames(totalLRpval)<-c("time","Lloc","Lgen","Lexp","Rloc","Rgen","Rexp","type","pval")
  lexpsTinp<-lexpsT[1:nrow(lexps)] 
  colnames(lexpsT)[(nrow(lexps)+1):length(colnames(lexpsT))]=c(colnames(ytyp),colnames(ytim))[c(colnames(ytyp),colnames(ytim)) !="id"]
  rexpsTinp<-rexpsT[1:nrow(rexps)]
  colnames(rexpsT)[(nrow(rexps)+1):length(colnames(rexpsT))]=c(colnames(ytyp),colnames(ytim))[c(colnames(ytyp),colnames(ytim)) !="id"]
  for(g in 1:length(totalLR$Lloc)){  
    lexpsotherinp<-as.data.frame(t(lexpsTinp[(as.character(lexpsT$Type)!=totalLR$Lloc[g] & as.character(lexpsT$Time)==totalLR$time[g]),colnames(lexpsTinp)==totalLR$Lgen[g]]))
    rexpsotherinp<-as.data.frame(t(rexpsTinp[(as.character(rexpsT$Type)!=totalLR$Rloc[g] & as.character(rexpsT$Time)==totalLR$time[g]),colnames(rexpsTinp)==totalLR$Rgen[g]]))
    lenl<-nrow(lexpsTinp[(as.character(lexpsT$Type)==totalLR$Lloc[g] & as.character(lexpsT$Time)==totalLR$time[g]),])
    lenr<-nrow(rexpsTinp[(as.character(rexpsT$Type)==totalLR$Rloc[g] & as.character(rexpsT$Time)==totalLR$time[g]),])
    
    pval<-select_base_pval(totalLR[g,],lexpsotherinp,rexpsotherinp,lenl,lenr) 
    totalLRpval[g,1:8]<-totalLR[g,1:8]
    totalLRpval[g,9]<-pval
  }
  kepLRpval<-na.omit(totalLRpval[totalLRpval$pval<=0.05,])
  write.csv(kepLRpval,paste(resloc,"singlenet.csv",sep=""),quote = F,row.names = F)
  
  for(t in 1:length(unique(kepLRpval$time))){
    kepLRpvalinp<-kepLRpval[kepLRpval$time==unique(kepLRpval$time)[t],]
    write.csv(kepLRpvalinp,paste(resloc,unique(kepLRpval$time)[t],"-singlenet.csv",sep=""),quote = F,row.names = F)
    single_network_plot(kepLRpvalinp,unique(kepLRpval$time)[t],usedcolo)
  }
}

developcorrelation<-function(ytim,usedcolo,resloc){ 
  data1<-read.csv(paste(resloc,"tmp/","data1.csv",sep=""))
  data2<-read.csv(paste(resloc,"tmp/","data2.csv",sep=""))
  sigres<-stat_LRcor(data1,data2,ytim)
  plot_LR_exps_label(data1,data2,sigres,ytim,usedcolo) 
}

comparednet<-function(controltime,treattime,usedcolo,resloc){
  kepLRpval<-read.csv(paste(resloc,"singlenet.csv",sep=""))
  if(treattime!="total"){
  kepLRpval1<-kepLRpval[kepLRpval$time %in% c(controltime,treattime),]
  }else{
  kepLRpval1<-kepLRpval  
  }
  compare_network_plot(kepLRpval1,controltime,usedcolo)
}

sankey<-function(time,type,usedcolo,mainlab,resloc){
  mergdat<-as.data.frame(prop.table(x = table(time,type), margin = 1))
  colnames(mergdat)<-c("Time","Type","Freq")
  mergdata<-mergdat[order(mergdat$Time),]
  s<-c()
  for(i in 1:length(unique(mergdata$Time))){
    if(i==1){
      s<-c(1:dim(mergdata[mergdata$Time==unique(mergdata$Time)[i],])[1])
    }else{
      s<-c(s,1:dim(mergdata[mergdata$Time==unique(mergdata$Time)[i],])[1])
    }
  }
  mergdata$subject<-s
  library(stringr)
  mergdata$Time<-factor(mergdata$Time,levels=str_sort(unique(mergdata$Time), numeric = TRUE))
  mergdata$Type<-factor(mergdata$Type,levels =unique(mergdata$Type))
  col1<-usedcolo[levels(mergdata$Type)]
  p<-ggplot(mergdata,
         aes(x = Time, stratum = Type, alluvium = subject,
             y = Freq,
             fill = Type)) + 
    scale_x_discrete(expand = c(.1, .1)) +
    geom_flow() +
    geom_stratum(alpha = .5) +
    scale_fill_manual(values = col1)+
    ggtitle(mainlab)+theme_classic()
  pdf(paste(resloc,"sankey.pdf",sep=""))
  print(p+theme(axis.text.x =element_text(angle=30)))
  dev.off()
}

definecolo<-function(tsccnetcol,type,resloc,mainlab){
  usedcolo<-tsccnetcol[1:length(unique(type))]
  names(usedcolo)<-unique(type)
  save(usedcolo,file=paste(resloc,mainlab,"-usedcolo.rdata",sep=""))
  return(usedcolo)
}

importance_plot<-function(resloc,ifscale,feature,columnorder){
  files<-list.files(path = resloc,pattern = "^Compared.*detail.csv")
  filenam<-sub("^Compared.*in ","",sub("netdetail.csv","",files))
  for(i in 1:length(files)){
    if(i==1){
      dat<-read.csv(paste(resloc,files[i],sep=""),row.names = 1)
      if(ifscale=="Column" |ifscale=="column" ){
        datas<-as.data.frame(scale(dat))
      }else{
        datas<-dat 
      }
    }else{
      dat0<-read.csv(paste(resloc,files[i],sep=""),row.names = 1)
      if(ifscale=="Column" | ifscale=="column"){
        datas0<-as.data.frame(scale(dat0))
      }else{
        datas0<-dat0 
      }
      datas<-rbind(datas,datas0)
    } 
  }
  print(paste("Have features: ",paste(colnames(datas),collapse = ", "),sep=""))
  datas$Celltype<-sub(".*:","",rownames(datas))
  datas$Time<-sub("^Compared.*in ","",sub("netdetail.csv","",sub(":.*","",rownames(datas))))

  inp0<-reshape2::dcast(datas,Celltype~Time,value.var=feature)
  inp<-inp0[,-1]
  rownames(inp)<-inp0[,1]
  
  sdnum<-apply(inp, 1, sd)
  inp1<-inp[sdnum !=0 & !is.na(sdnum),]
  if(columnorder=="None"){
  inps<-inp1  
  }else{
  inps<-inp1[,columnorder]
  }
  pdf(paste(resloc,"Importance-",feature,"-",ifscale,"scaled",".pdf",sep=""))
  if(ifscale=="row" | ifscale=="Row"){
    p<-pheatmap(inps,scale = "row",na_col = "gray",cluster_cols = F,treeheight_row = 0,angle_col=45, color = colorRampPalette(c("navy", "white", "firebrick3"))(100),border_color = NA,main = feature) 
    print(p)
  }else{
    p<-pheatmap(inps,na_col = "gray",cluster_cols = F,treeheight_row = 0,angle_col=45, color = colorRampPalette(c("navy", "white", "firebrick3"))(100),border_color = NA,main = feature) 
    print(p)
  }
  dev.off()
}

