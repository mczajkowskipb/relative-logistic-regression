# Frozen Stage H1 fix1: unambiguous budget and duplicate guard.
# All tests use only base R, and never read outcomes or V/E.
critH1_need_pool <- function(pool,method,budget=320L) {
  allowed <- if (identical(method,"EXHAUSTIVE32")) 496L else as.integer(budget)
  valid_shape <- is.matrix(pool) && nrow(pool)==2L &&
     is.numeric(pool) && all(is.finite(pool))
  if(!isTRUE(valid_shape))
    stop("STAGEH1_STOP: candidate pool matrix shape invalid",call.=FALSE)
  K <- ncol(pool)
  if(K>allowed || K<1L)
    stop(sprintf(
      "STAGEH1_STOP: pool budget exceeded method=%s K=%d allowed=%d",
      method,K,allowed),call.=FALSE)
  if(any(pool[1L,]>=pool[2L,]))
    stop("STAGEH1_STOP: noncanonical unordered candidate pairs",call.=FALSE)
  if(anyDuplicated(paste(pool[1L,],pool[2L,],sep=":"))>0L)
    stop("STAGEH1_STOP: duplicated unordered candidate pair",call.=FALSE)
  invisible(TRUE)
}
