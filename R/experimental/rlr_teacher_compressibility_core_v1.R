# Exact bounded rank-tree RLR compressibility DIAGNOSTIC, not a new method.
# Input p_oof: honest D-only crossfit original RLR20 probabilities.
source("R/experimental/rlr_global_logit_variance_tree_core_v5.R")
rc_need <- rrv5_need
rc_fit <- function(D,y,pairs,p_oof,cap=3L,min_leaf=4L,goal=c("TEACHER","LABEL")){
 goal<-match.arg(goal)
 rc_need(is.factor(y)&&nlevels(y)==2L&&nrow(D)==length(y),"binary D labels")
 rc_need(cap %in% 0:3 && min_leaf>=2,"invalid fixed tree budget")
 rc_need(length(p_oof)==nrow(D)&&all(is.finite(p_oof))&&
           all(p_oof>0&p_oof<1),"teacher must be valid crossfitted D OOF")
 X<-rrv5_matrix(D,pairs);n<-nrow(X)
 yy<-as.integer(y==levels(y)[2]);logit<-qlogis(pmin(.99999,pmax(.00001,p_oof)))
 prior<-(sum(yy)+1)/(n+2)
 py<-function(ix)(sum(yy[ix])+4*prior)/(length(ix)+4)
 pt<-function(ix)plogis(mean(logit[ix]))
 loss<-function(ix){
   if(goal=="TEACHER")return(sum((logit[ix]-mean(logit[ix]))^2))
   p<-py(ix);-sum(yy[ix]*log(p)+(1-yy[ix])*log1p(-p))
 }
 root<-seq_len(n)
 base<-loss(root)
 subbest<-function(ix,excluded){
   out<-list(feature=NA_integer_,gain=0)
   if(length(ix)<2*min_leaf)return(out)
   for(j in seq_len(ncol(X))){
     if(j==excluded)next
     a<-ix[X[ix,j]==0L];b<-ix[X[ix,j]==1L]
     if(min(length(a),length(b))<min_leaf)next
     g<-loss(ix)-loss(a)-loss(b)
     if(g>out$gain+1e-12)out<-list(feature=j,gain=g,a=a,b=b)
   }
   out
 }
 win<-list(feature=NA_integer_,cost=base)
 if(cap>0L)for(j in seq_len(ncol(X))){
   a<-which(X[,j]==0L);b<-which(X[,j]==1L)
   if(min(length(a),length(b))<min_leaf)next
   c0<-c1<-list(feature=NA_integer_,gain=0)
   if(cap>=2L){
     ca<-subbest(a,j);cb<-subbest(b,j)
     if(cap==2L){
       if(ca$gain>=cb$gain && ca$gain>0)c0<-ca
       if(cb$gain>ca$gain && cb$gain>0)c1<-cb
     }else{c0<-ca;c1<-cb}
   }
   val<-loss(a)+loss(b)-c0$gain-c1$gain
   if(val<win$cost-1e-12)win<-list(feature=j,cost=val,a=a,b=b,c0=c0,c1=c1)
 }
 nodes<-list()
 add<-function(ix,dep,plan=NULL){
   id<-length(nodes)+1L
   f<-if(is.null(plan))NA_integer_ else plan$feature
   nodes[[id]]<<-list(depth=dep,n=length(ix),f=f,kids=integer(),
                     prob_label=py(ix),prob_teacher=pt(ix))
   if(!is.na(f)){
     l<-add(plan$a,dep+1L,if(dep==0L)plan$c0 else NULL)
     r<-add(plan$b,dep+1L,if(dep==0L)plan$c1 else NULL)
     nodes[[id]]$kids<<-c(l,r)
   }
   id
 }
 add(root,0L,win)
 used<-sort(unique(vapply(nodes,function(z)z$f,integer(1))),na.last=NA)
 slim<-pairs[used,c("gene_left","gene_right"),drop=FALSE];rownames(slim)<-NULL
 for(i in seq_along(nodes))if(!is.na(nodes[[i]]$f))
   nodes[[i]]$f<-match(nodes[[i]]$f,used)
 internal<-sum(vapply(nodes,function(z)!is.na(z$f),logical(1)))
 rc_need(internal<=cap&&internal<=3L&&length(nodes)==2L*internal+1L,
         "invalid finite optimizer tree")
 structure(list(nodes=nodes,pairs=slim,goal=goal,cap=cap,
   internal=internal,leaves=internal+1L,
   n_unique_pairs=length(used),
   n_unique_genes=length(unique(c(slim$gene_left,slim$gene_right))),
   train_loss=win$cost,root_loss=base,
   train_fraction=if(base>1e-12)win$cost/base else 0,
   n_D=n),class="rlr_compressibility_diagnostic")
}
rc_predict<-function(fit,X){
 R<-rrv5_ranks(X);P<-fit$pairs
 rc_need(all(c(P$gene_left,P$gene_right)%in%colnames(R)),
         "missing trained gene universe")
 n<-nrow(R);pl<-pt<-numeric(n);count<-integer(n)
 path<-character(n);leaf<-integer(n)
 for(i in seq_len(n)){
   id<-1L;h<-character()
   repeat{
     node<-fit$nodes[[id]]
     if(is.na(node$f)){
       pl[i]<-node$prob_label;pt[i]<-node$prob_teacher;leaf[i]<-id;break
     }
     f<-node$f
     branch<-as.integer(R[i,P$gene_right[f]]>R[i,P$gene_left[f]])+1L
     count[i]<-count[i]+1L
     h<-c(h,paste0(P$gene_left[f],">",P$gene_right[f],":",branch))
     id<-node$kids[branch]
   }
   path[i]<-paste(h,collapse=" / ")
 }
 rc_need(all(count<=2L)&&all(pl>0&pl<1)&&all(pt>0&pt<1),
         "invalid autonomous student prediction")
 list(prob_label=pl,prob_teacher=pt,checks=count,path=path,leaf=leaf)
}
