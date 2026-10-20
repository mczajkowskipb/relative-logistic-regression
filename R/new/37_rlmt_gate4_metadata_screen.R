# GATE 4: blind metadata/integrity screen for the six uninspected RLMT
# CompCancer datasets. NO algorithm fitting or predictive metrics.
# Run: Rscript --vanilla R/new/37_rlmt_gate4_metadata_screen.R
# Eligibility fixed before inspection: docs/RLMT_GATE4_METADATA_RULES.md

if (!requireNamespace("foreign", quietly = TRUE))
  stop("foreign package is required to inspect ARFF source")
source_dir <- "/opt/2026/RLR/data/CompCancer_truth_9x10"
manifest <- "checksums/DATA_SHA256_EXPECTED.txt"
output <- "results/rlmt_gate4_metadata_screen_v1"

datasets <- c(
  "alizadeh-2000-v1", "armstrong-2002-v1", "bittner-2000",
  "chen-2002", "golub-1999-v1", "gordon-2002",
  "pomeroy-2002-v1", "shipp-2002-v1", "singh-2002"
)
development <- c("alizadeh-2000-v1", "chen-2002", "singh-2002")
remaining <- setdiff(datasets, development)
fail <- function(...) stop(..., call. = FALSE)
if (dir.exists(output)) fail("STOP: Refusing to overwrite ", output)
if (!dir.exists(source_dir)) fail("Missing CompCancer source directory")
if (!file.exists(manifest)) fail("Missing canonical SHA256 manifest")
if (!nzchar(Sys.which("sha256sum"))) fail("sha256sum executable not found")

files <- list.files(source_dir, pattern = "[.]arff$", full.names = FALSE)
expected_file_names <- paste0(datasets, ".arff")
if (!setequal(files, expected_file_names) ||
    length(files) != length(expected_file_names)) {
  fail("STOP: CompCancer inventory does not match fixed nine-file scope: ",
       paste(files, collapse = ", "))
}

lines <- readLines(manifest, warn=FALSE)
expected_sha <- character()
for (line in lines) {
  if (!grepl("data/CompCancer/[^[:space:]]+[.]arff$", line)) next
  parts <- strsplit(trimws(line), "[[:space:]]+")[[1L]]
  if (length(parts)!=2L || !grepl("^[[:xdigit:]]{64}$",parts[1L]))
    fail("Unrecognized canonical checksum manifest entry: ",line)
  basename_file <- basename(parts[2L])
  if (basename_file %in% names(expected_sha))
    fail("Duplicate manifest filename: ", basename_file)
  expected_sha[basename_file] <- tolower(parts[1L])
}
if (!all(expected_file_names %in% names(expected_sha)))
  fail("Missing canonical checksums for some of the nine data files")

sha256 <- function(path) {
  out <- suppressWarnings(system2(
    "sha256sum", args = shQuote(normalizePath(path)),
    stdout = TRUE, stderr = TRUE))
  if (length(out)!=1L ||
      !grepl("^[[:xdigit:]]{64}[[:space:]]+",out)) {
    fail("Failed to compute SHA256 for ",path)
  }
  tolower(sub("^([[:xdigit:]]{64}).*$","\\1",out))
}

audit_one <- function(dataset) {
  file_name <- paste0(dataset,".arff")
  p <- file.path(source_dir,file_name)
  actual_sha <- sha256(p)
  check_sha <- identical(actual_sha,unname(expected_sha[[file_name]]))
  data <- tryCatch(foreign::read.arff(p),error=function(e)e)
  if (inherits(data,"error")) {
    return(data.frame(
      dataset=dataset, development_used=dataset %in% development,
      file_bytes=as.numeric(file.info(p)$size), sha256=actual_sha,
      sha256_matches=check_sha, parse_ok=FALSE,
      parse_error=conditionMessage(data), class_col=NA_character_,
      class_levels=NA_character_, class_counts=NA_character_,
      n_samples=NA_integer_, min_class_n=NA_integer_,
      n_numeric_features=NA_integer_, n_other_features=NA_integer_,
      nonfinite_numeric_cells=NA_integer_,
      duplicate_numeric_rows=NA_integer_,
      patient_id_verified=FALSE, candidate_eligible=FALSE,
      eligibility_reason="ARFF_PARSE_FAILED",stringsAsFactors=FALSE))
  }
  candidates <- c("class","Class","CLASS","label","Label","y","Y")
  hits <- intersect(candidates,colnames(data))
  identified <- length(hits)==1L
  class_name <- if (identified) hits[[1L]] else NA_character_
  raw_y <- if (identified) data[[class_name]] else NULL
  y <- if (identified) factor(raw_y) else NULL
  class_ok <- identified && !anyNA(y) && length(levels(y))==2L
  class_n <- if (class_ok) as.integer(table(y)) else integer()
  other <- if (identified) setdiff(names(data),class_name) else names(data)
  features <- data[,other,drop=FALSE]
  numeric <- vapply(features,is.numeric,logical(1))
  X <- features[,numeric,drop=FALSE]
  n_nonfinite <- if (ncol(X)) {
    as.integer(sum(vapply(X,function(column)
      sum(!is.finite(column)),integer(1))))
  } else 0L
  # Complete-row duplication uses all numeric model predictors, independent
  # of labels. No data is removed; report only.
  n_duplicates <- if (ncol(X)) {
    as.integer(sum(duplicated(X)))
  } else NA_integer_
  # Existing RLR read_xy_file() synthesizes row names for ARFF without
  # explicit identifiers. ARFF contents here cannot establish patient IDs.
  eligible <- !(dataset %in% development) &&
    check_sha && class_ok && nrow(data)>=80L &&
    min(class_n)>=25L && ncol(X)>=20L &&
    n_nonfinite==0L && !is.na(n_duplicates) && n_duplicates==0L
  reasons <- character()
  if (dataset %in% development) reasons <- c(reasons,"USED_FOR_RLMT_DEVELOPMENT")
  if (!check_sha) reasons <- c(reasons,"SHA256_MISMATCH")
  if (!identified) reasons <- c(reasons,"CLASS_COLUMN_ABSENT_OR_AMBIGUOUS")
  if (identified && !class_ok) reasons <- c(reasons,"NON_BINARY_OR_MISSING_CLASS")
  if (nrow(data)<80L) reasons <- c(reasons,"N_LT_80")
  if (class_ok && min(class_n)<25L) reasons <- c(reasons,"MIN_CLASS_LT_25")
  if (ncol(X)<20L) reasons <- c(reasons,"NUMERIC_FEATURES_LT_20")
  if (n_nonfinite>0L) reasons <- c(reasons,"NONFINITE_NUMERIC_VALUES")
  if (!is.na(n_duplicates) && n_duplicates>0L)
    reasons <- c(reasons,"DUPLICATED_NUMERIC_ROWS")
  if (eligible) reasons <- "METADATA_SCREEN_ELIGIBLE"
  data.frame(
    dataset=dataset, development_used=dataset %in% development,
    file_bytes=as.numeric(file.info(p)$size), sha256=actual_sha,
    sha256_matches=check_sha, parse_ok=TRUE, parse_error="",
    class_col=class_name,
    class_levels=if (identified) paste(levels(y),collapse=" | ") else NA_character_,
    class_counts=if (identified) paste(as.integer(table(y)),collapse=" | ") else NA_character_,
    n_samples=nrow(data),
    min_class_n=if(class_ok)min(class_n)else NA_integer_,
    n_numeric_features=ncol(X), n_other_features=sum(!numeric),
    nonfinite_numeric_cells=n_nonfinite,
    duplicate_numeric_rows=n_duplicates,
    patient_id_verified=FALSE, candidate_eligible=eligible,
    eligibility_reason=paste(reasons,collapse=";"),
    stringsAsFactors=FALSE)
}
cat("=== GATE 4 METADATA SCREEN (NO MODELS) ===\n")
rows <- lapply(datasets,function(ds) {
  cat("READ ",ds,"\n",sep="")
  audit_one(ds)
})
tab <- do.call(rbind,rows)
if (nrow(tab)!=9L || anyDuplicated(tab$dataset))
  fail("Bad inventory output")
if (!all(tab$sha256_matches) || !all(tab$parse_ok)) {
  cat("INTEGRITY_OR_PARSE_FAILURE=TRUE\n")
}

cat("\n=== METADATA SUMMARY ===\n")
print(tab[,c("dataset","development_used","sha256_matches",
             "n_samples","class_levels","class_counts",
             "n_numeric_features","nonfinite_numeric_cells",
             "duplicate_numeric_rows","candidate_eligible",
             "eligibility_reason")],row.names=FALSE,right=FALSE)

quality_ok <- all(tab$sha256_matches) && all(tab$parse_ok)
eligible_set <- tab$dataset[tab$candidate_eligible]
cat("\nMETADATA_SCREEN_ALL_HASHES=",if(quality_ok)"PASS" else "FAIL","\n",sep="")
cat("UNINSPECTED_CANDIDATES=",length(remaining),"\n",sep="")
cat("METADATA_ELIGIBLE_COUNT=",length(eligible_set),"\n",sep="")
cat("METADATA_ELIGIBLE_SET=",if(length(eligible_set))
  paste(eligible_set,collapse=";") else "NONE","\n",sep="")
cat("PATIENT_LEVEL_INDEPENDENCE=UNKNOWN\n")
cat("EXTERNAL_VALIDATION=NO_HISTORICAL_COMPCANCER_ONLY\n")

# Writes summary only: source ARFF, existing experiment results untouched.
dir.create(output,recursive=TRUE,showWarnings=FALSE)
write.csv(tab,file.path(output,"metadata_inventory.csv"),row.names=FALSE)
writeLines(c(
  "GATE4 metadata screen. No training, CV, predictions or test-label selection.",
  "SHA256 checked against the locked RLR publication canonical CompCancer manifest.",
  paste("Eligibility: all six untouched RLMT-development datasets satisfying",
        "binary classes; n>=80; minority>=25; numeric p>=20;",
        "zero nonfinite numeric cells; zero exact numeric row duplicates;",
        "canonical checksum match."),
  paste("Integrity:", if (quality_ok) "PASS" else "FAIL"),
  paste("Metadata-eligible:",if (length(eligible_set))
        paste(eligible_set,collapse=",") else "NONE"),
  "No patient IDs verified in processed CompCancer matrices.",
  "Not a newly sourced external cohort; cannot claim prospective external validation.",
  "Any later model evaluation is a separate explicitly frozen run."
),file.path(output,"README.txt"))

if (!quality_ok) fail("GATE4_METADATA_SCREEN=FAIL_INTEGRITY")
cat("GATE4_METADATA_SCREEN=PASS\n")
cat("REPORT=",output,"\n",sep="")
