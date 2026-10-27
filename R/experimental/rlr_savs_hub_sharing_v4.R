# SAVS4: root/leaf gene hub SHARING (not same pair duplication), motivated
# by user's x1 > {x2,x3,x4} hypothesis. Root share is allowed; child rules
# remain within-branch gene-disjoint. All candidate ranking is samplewise.
source("R/experimental/rlr_savs_honest_lookahead_v3.R")
sv4_need<-rg_need

sv4_exclude_root_exact_pairs<-function(pool,root){
 if(is.null(root))return(pool)
 key<-function(a,b)paste(pmin(a,b),pmax(a,b),sep="|")
 prior<-key(rep(root$target,length(root$anchors)),root$anchors)
 curr<-key(as.character(pool$gene_left),as.character(pool$gene_right))
 pool[!curr %in% prior,,drop=FALSE]
}
sv4_fit <- function(D,y,pool,root,mode="RLR10"){
 sv4_need(!is.null(root),"hub-sharing requires a genuine root")
 R<-rg_ranks(D);g<-os_gate(R,root)
 p<-(sum(sv_truth(y))+2)/(length(y)+4)
 m<-sv_m(ncol(D),mode)
 panel_pool<-sv4_exclude_root_exact_pairs(pool,root)
 chosen<-vector("list",2L);fits<-vector("list",2L)
 for(b in 0:1){
  ids<-which(g==b)
  # Crucial difference from SAVS1: DO NOT ban x_root or root anchors
  # from other, NEW child pairs. The exact root pair is still banned.
  chosen[[b+1L]]<-sv_choose_rules(R,y,panel_pool,ids,m,
                       max_per_class=2L,excluded=character())
  fits[[b+1L]]<-sv_bucket_fit(R,y,ids,chosen[[b+1L]],
                              m,global_prior=p)
 }
 structure(list(root=root,rule_panels=chosen,fits=fits,
    margin=m,prior=p,
    training_branches=c(sum(g==0L),sum(g==1L)),
    n_rules=vapply(chosen,nrow,integer(1)),
    n_path_gene_overlap=vapply(chosen,function(tab){
      if(!nrow(tab))return(0L)
      length(intersect(unique(c(tab$gene_left,tab$gene_right)),
               c(root$target,root$anchors)))
    },integer(1))),
    class="SAVS4")
}
sv4_choose<-function(D,y,mode="RLR10",max_roots=12L){
 ids<-sv3_split(y)
 S<-D[ids$S,,drop=FALSE];H<-D[ids$H,,drop=FALSE]
 ys<-factor(as.character(y[ids$S]),levels=levels(y))
 yh<-factor(as.character(y[ids$H]),levels=levels(y))
 fallback<-list(accepted=FALSE,root=NULL,mode=mode,
     n_S=nrow(S),n_H=nrow(H),proposed=0L,
     H_globalF1=NA_real_,H_treeF1=NA_real_,
     H_globalBA=NA_real_,H_treeBA=NA_real_,
     reason="NO_SUPPORT")
 if(min(table(ys))<2L || min(table(yh))<1L)return(fallback)
 poolS<-compute_top_pairs_switchbox(S,ys,top_n=1000L,buffer_mult=5L)
 options<-sv2_roots(S,ys,poolS,max_hubs=8L,max_roots=max_roots)
 fallback$proposed<-length(options)
 control<-sv_global4(S,ys,poolS,sv_m(ncol(S),mode))
 baseline<-sv2_metrics(yh,sv_predict(control,H)$p)
 fallback$H_globalF1<-baseline$macroF1
 fallback$H_globalBA<-baseline$BA
 if(!length(options))return(fallback)
 tested<-lapply(options,function(c){
   model<-sv4_fit(S,ys,poolS,c$root,mode)
   z<-sv2_metrics(yh,sv_predict(model,H)$p)
   list(root=c$root,positiveF1=z$positiveF1,
        macroF1=z$macroF1,BA=z$BA,LL=z$LL)
 })
 order_idx<-order(-vapply(tested,function(x)x$macroF1,numeric(1)),
                  -vapply(tested,function(x)x$BA,numeric(1)),
                   vapply(tested,function(x)x$LL,numeric(1)),
                   vapply(tested,function(x)length(x$root$anchors),integer(1)),
                   seq_along(tested))
 win<-tested[[order_idx[1L]]]
 fallback$H_treeF1<-win$macroF1
 fallback$H_treeBA<-win$BA
 fallback$accepted<-isTRUE(win$macroF1>=baseline$macroF1+.01 &&
                            win$BA>=baseline$BA-.02)
 fallback$root<-win$root
 fallback$reason<-if(fallback$accepted)
  "S_H_HUB_SHARE_ACCEPTED" else "S_H_HUB_SHARE_GLOBAL_FALLBACK"
 fallback
}
sv4_refit<-function(D,y,pool,selection){
 if(isTRUE(selection$accepted))
   sv4_fit(D,y,pool,selection$root,selection$mode)
 else
   sv_global4(D,y,pool,sv_m(ncol(D),selection$mode))
}
