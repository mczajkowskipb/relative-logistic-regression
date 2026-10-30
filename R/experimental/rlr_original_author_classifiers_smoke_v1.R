# Unmodified OFFICIAL models on authors' original example matrices.
# This is an inference/reference-integrity check, NOT an RLR performance benchmark.
args<-commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
root<-normalizePath(args[1],mustWork=TRUE)
out<-file.path("results","rlr-original-author-classifier-smoke-v1")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
stopifnot(file.exists(file.path(root,"PurIST","functions.R")),
          file.exists(file.path(root,"PurIST","fitteds_public_2019-02-12.Rdata")),
          file.exists(file.path(root,"PurIST","example.Rdata")),
          file.exists(file.path(root,"DeCAF","decaf_functions.R")),
          file.exists(file.path(root,"DeCAF","decaf_classifier.Rdata")),
          file.exists(file.path(root,"DeCAF","example.rds")))
p_env<-new.env(parent=globalenv())
load(file.path(root,"PurIST","fitteds_public_2019-02-12.Rdata"),envir=p_env)
load(file.path(root,"PurIST","example.Rdata"),envir=p_env)
sys.source(file.path(root,"PurIST","functions.R"),envir=p_env)
stopifnot(exists("classifs",p_env),exists("dat",p_env),
          exists("apply_classifier",p_env))
pm<-p_env$classifs[[1L]]
pdat<-p_env$dat
pr<-p_env$apply_classifier(data=pdat,classifier=pm)
stopifnot(is.data.frame(pr), "Pred_prob_basal" %in% names(pr),
          nrow(pr)==ncol(pdat),
          all(is.finite(pr$Pred_prob_basal)),
          all(pr$Pred_prob_basal>=0 & pr$Pred_prob_basal<=1))
write.csv(data.frame(model="PurIST_original_frozen",sample_id=colnames(pdat),
                     predicted_probability=as.numeric(pr$Pred_prob_basal),
                     outcome_name="basal-like PDAC subtype"),
          file.path(out,"purist_official_example_predictions.csv"),row.names=FALSE)
cat(sprintf("ORIGINAL_PurIST_FROZEN_AUTHOR_CODE_EXAMPLE=PASS samples=%d pairs=%d\n",
            ncol(pdat),nrow(pm$TSPs)))
d_env<-new.env(parent=globalenv())
load(file.path(root,"DeCAF","decaf_classifier.Rdata"),envir=d_env)
sys.source(file.path(root,"DeCAF","decaf_functions.R"),envir=d_env)
ddat<-readRDS(file.path(root,"DeCAF","example.rds"))
stopifnot(exists("decaf_classifier",d_env),exists("apply_decaf",d_env))
dm<-d_env$decaf_classifier$classifier2
dr<-d_env$apply_decaf(data=ddat,classifier=dm)
stopifnot(is.data.frame(dr),"DeCAF_prob" %in% names(dr),
          nrow(dr)==ncol(ddat),
          all(is.finite(dr$DeCAF_prob)),
          all(dr$DeCAF_prob>=0 & dr$DeCAF_prob<=1))
write.csv(data.frame(model="DeCAF_original_frozen",sample_id=colnames(ddat),
                     predicted_probability=as.numeric(dr$DeCAF_prob),
                     outcome_name="proCAF subtype"),
          file.path(out,"decaf_official_example_predictions.csv"),row.names=FALSE)
cat(sprintf("ORIGINAL_DeCAF_FROZEN_AUTHOR_CODE_EXAMPLE=PASS samples=%d pairs=%d\n",
            ncol(ddat),nrow(dm$TSPs)))
cat("ORIGINAL_PurIST_DeCAF_REFERENCE_SMOKE_NOT_A_HEADTOHEAD=PASS\n")
