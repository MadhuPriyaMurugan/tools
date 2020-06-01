#!/ usr/bin/Rscript

###
### Pairwise BUB and Jaccard index for pairwise mouse samples - tissue version
###

###
### Using nucleotide sequences to determine similarity.
###

###
### NOTES
###
### THIS IS ONLY FOR WITHIN BATCH
### THIS IS ONLY TO COMPARE DIFFERENT TISSUE SAMPLES FROM A SINGLE MOUSE

### Dependencies
suppressMessages(library(data.table))
install.packages("optparse")
library(optparse)
suppressMessages(library(VennDiagram))
library(MASS)
library(scales)
suppressMessages(library(xlsx))
source("~/my_tool_repos/WesPersonal/utilityFxns.R")

####################
### COMMAND LINE ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####################

### Make list of options
optlist <- list(
  make_option(
    c("-i", "--inputFile"),
    type = "character",
    help = "full_clones.txt file created by groupClones.R"
  ),
  make_option(
    c("-m", "--meta"),
    type = "character",
    help = "Metadata file containing treatment designations at a minimum. Require: Sample | Tissue | Mouse/Animal | Treatment (optional)"
  ),
  make_option(
    c("-p", "--pairs"),
    type = "character",
    help = "Metadata column name used to check correct pairing. Could be sample number if separate batches, or mouse/animal ID if same batch. Comma-sep, no spaces if more than one."
  ),
  make_option(
    c("-c", "--compare"),
    type = "character",
    default = "tissue",
    help = "Metadata column used to compare. 'tissue' (default) will compare by tissue. 'treatment' will compare by treatment."
  ),
  make_option(
    c("-s", "--subset"),
    type = "character",
    help = "Which samples should be extracted. Should either be in format 'm:n' or 'l,m,n,o,p' or 'm:n,o:p'."
  ),
  make_option(
    c("-f", "--freqGroup"),
    type = "character",
    help = "Which frequency groups should be included in the comparison. Two comma-separated, no space vectors that are themselves
    separated by a semicolon. First element corresponds to the first value in 'compare' column of metadata. The following example
    will compare hyperexpanded to small, medium, and large: 'hyperexpanded;small,medium,large'. Can also use 'all' for all groups."
  ),
  make_option(
    c("-o", "--outDir"),
    type = "character",
    help = "file path to directory for writing output files"
  ),
  make_option(
    c("-n", "--wkbkName"),
    type = "character",
    help = "Name of output excel workbook."
  ),
  make_option(
    c("-S", "--sheetName"),
    type = "character",
    help = "Name of workbook sheet."
  ),
  make_option(
    c("-d", "--debug"),
    type = "logical",
    default = FALSE,
    help = "Logical. TRUE - print session info and extra output to help with debugging. Also do not write output (tables and images). FALSE - normal output and write output files (tables and images)."
  ),
  make_option(
    c("-l", "--log"),
    type = "logical",
    default = FALSE,
    help = "Logical. TRUE - output session info. FALSE - do not output session info. If debug is TRUE and log is FALSE, then session info will be printed to STDOUT. If neither are set, no output."
  )
)



### Parse commandline
p <- OptionParser(usage = "%prog -i inputDir -m meta -p pair -c compare -s subset -f freqGroup -o outDir -n wkbkName -S sheetName -d debug -l log",
                  option_list = optlist)
args <- parse_args(p)
opt <- args$options

### Commands
inputFile_v <- args$inputFile
metaFile_v <- args$meta
pair_v <- args$pair
compare_v <-args$compare
subset_v <- args$subset
freqGroup_v <- args$freqGroup
outDir_v <- args$outDir
wkbkName_v <- args$wkbkName
sheetName_v <- args$sheetName
debug_v <- args$debug
log_v <- args$log

source("~/my_tool_repos/WesPersonal/utilityFxns.R")

### LIB180920LC
# inputFile_v <- "~/OHSU/tcr_spike/data/LIB180920LC/freqGroups/LIB180920LC_full_clones.txt"
# metaFile_v <- "~/OHSU/tcr_spike/data/LIB180920LC/meta/meta.txt"
# pair_v <- "Mouse"
# compare_v <- "Tissue"
# outDir_v <- "~/OHSU/tcr_spike/data/LIB180920LC/jaccard/"
# wkbkName_v <- "LIB180920LC_bloodAll_TumorSmall_Jaccard.xlsx"
# debug_v <- T
# log_v <- T
# subset_v <- NULL
# freqGroup_v <- "All;Small"
# sheetName_v <- NULL

# ### LIB181105LC
# inputFile_v <- "~/OHSU/tcr_spike/data/LIB181105LC/freqGroups/LIB181105LC_full_clones.txt"
# metaFile_v <- "~/OHSU/tcr_spike/data/LIB181105LC/meta/meta.txt"
# pair_v <- "Mouse"
# compare_v <- "Tissue"
# outDir_v <- "~/OHSU/tcr_spike/data/LIB181105LC/jaccard/"
# wkbkName_v <- "LIB181105LC_lymphLungJaccard.xlsx"
# debug_v <- T
# log_v <- T
# subset_v <- NULL
# freqGroup_v <- "Hyperexpanded;All"
# sheetName_v <- NULL

### Handle subset argument
if (!is.null(subset_v)) {
  subset_v <- unlist(strsplit(subset_v, ","))
  subset_v <- unlist(sapply(subset_v, function(x) {
    y <- strsplit(x, ":")[[1]]
    z <- as.numeric(y[1]):as.numeric(y[2])
    return(z)
  }, simplify = F, USE.NAMES = F))
}

### Handle sheet argument
sheetName_v <- ifelse(is.null(sheetName_v), "Sheet1", sheetName_v)

### Handle freqGroup argument
if (!is.null(freqGroup_v)) {
  freqGroup_lsv <- as.list(strsplit(freqGroup_v, split = ";")[[1]])
  freqGroup_lsv <- lapply(freqGroup_lsv, function(x) {
    if (x %in% c("all", "All")) {
      y <- c("Rare", "Small", "Medium", "Large", "Hyperexpanded")
    } else {
      y <- unlist(strsplit(x, split = ','))
    }
    y
  })
}

### Check arguments
print("Wrangled subset:")
print(subset_v)
print("")
print("Current input file:")
print(inputFile_v)

#############
### SETUP ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#############

### Print commands
if (log_v){
  returnSessionInfo(args_lsv = args, out_dir_v = outDir_v)
} else {
  if (debug_v){
    returnSessionInfo(args_lsv = args)
  } # fi
} # fi

### Get data
meta_dt <- fread(metaFile_v)
inputData_dt <- fread(inputFile_v)

### Get column names
sampleCol_v <- grep("[Ss]ample", colnames(meta_dt), value = T)
animalCol_v <- grep("[Mm]ouse|[Aa]nimal", colnames(meta_dt), value = T)
treatCol_v <- grep("[Tt]reatment", colnames(meta_dt), value = T)
tissueCol_v <- grep("[Tt]issue", colnames(meta_dt), value = T)

### Subset metadata
if (!is.null(subset_v)) {
  keepRows_v <- subset_v
  meta_dt <- meta_dt[keepRows_v,]
}

print("Current meta:")
print(meta_dt)

## Get pair column name
## This is the column that identifies which mice are pairs
pairCol_v <- unlist(sapply(pair_v, function(x){
  y <- grep(x, colnames(meta_dt), value = T)
  if (length(y) == 0) {
    y <- grep(simpleCap(pair_v), colnames(meta_dt), value = T)
  }
  return(y)
}, simplify = F))

### Get different tissues that will be compared
tissues_v <- unique(meta_dt[,get(tissueCol_v)])

### Name freq group divisions with tissues
names(freqGroup_lsv) <- tissues_v

### Get clone identity columns
seqCol_v <- "nSeqCDR3"
vCol_v <- grep("V Segments|V segments", colnames(inputData_dt), value = T)
jCol_v <- grep("J Segments|J segments", colnames(inputData_dt), value = T)

### Add compare column
inputData_dt$compare <- paste(inputData_dt[[vCol_v]], inputData_dt[[seqCol_v]], inputData_dt[[jCol_v]], sep = "_")

###################
### CHECK INPUT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###################

### Determine how many replicates each tissue will have
byMouse_dt <- meta_dt[,.N, by = pairCol_v]
byTissue_dt <- meta_dt[,.N, by = tissueCol_v]

######################
### PREPARE OUTPUT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
######################

## Make data.tables
output_dt <- as.data.table(matrix(nrow = nrow(meta_dt), ncol = 8))
colnames(output_dt) <- c("PairSamples", "BUB_Value", "J_Value", "TotalClones",
                         "ClonesA", "ClonesB", "Shared", "Union")

## Change to character or numeric
charCols_v <- colnames(output_dt)[1]
numCols_v <- colnames(output_dt)[2:ncol(output_dt)]
for (col_v in charCols_v) set(output_dt, j = col_v, value = as.character(output_dt[[col_v]]))
for (col_v in numCols_v) set(output_dt, j = col_v, value = as.numeric(output_dt[[col_v]]))

## Add meta info
output_dt <- cbind(meta_dt, output_dt)

############################
### PERFORM CALCULATIONS ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
############################

### Each row in the `totalPairs_dt` data.table corresponds to a unique mouse tissue pair
### Iterate over those, using the experiment and mouse IDs to grab the appropriate batches
### and samples from the metadata. Then read them in and perform the calculations.

for (i in 1:byMouse_dt[,.N]){
  print(sprintf("Index %d", i))
  ## Get information
  currAnimal_v <- byMouse_dt[i,get(animalCol_v)]
  
  ## Subset metadata
  currMeta_dt <- meta_dt[get(animalCol_v) == currAnimal_v,]
  
  ## Read in data
  currData_lsdt <- list()
  for (j in 1:nrow(currMeta_dt)){
    ## Get sample and batch
    currTissue_v <- currMeta_dt[j,get(tissueCol_v)]
    currSamp_v <- currMeta_dt[j,get(sampleCol_v)]
    currName_v <- paste(currTissue_v, currSamp_v, currAnimal_v, sep = "_")
    ## Subset for specific sample
    currData_dt <- inputData_dt[get(sampleCol_v) == paste0("S", gsub("^S", "", currSamp_v)),]
    ## Subset for specific freq groups
    currGroups_v <- freqGroup_lsv[[currTissue_v]]
    currData_dt <- currData_dt[Div %in% currGroups_v,]
    ## Add file
    currData_lsdt[[currName_v]] <- currData_dt
  } # for j
  
  ## Get total number of clones in all samples
  totalClones_v <- unique(c(unlist(sapply(currData_lsdt, function(x) x$compare))))
  
  ## Get current sample
  currFirst_v <- names(currData_lsdt)[1]
  currFSample_v <- strsplit(currFirst_v, "_")[[1]][2]
  currFTissue_v <- strsplit(currFirst_v, "_")[[1]][1]
  currFMouse_v <- strsplit(currFirst_v, "_")[[1]][3]
  
  ## Get other sample
  currSecond_v <- names(currData_lsdt)[2]
  currSSample_v <- strsplit(currSecond_v, "_")[[1]][2]
  currSTissue_v <- strsplit(currSecond_v, "_")[[1]][1]
  currSMouse_v <- strsplit(currSecond_v, "_")[[1]][3]
  
  ## Get clones
  currFClones_v <- currData_lsdt[[currFirst_v]]$compare
  currSClones_v <- currData_lsdt[[currSecond_v]]$compare
  
  ## Get intersection, union, and absent
  currInt_v <- intersect(currFClones_v, currSClones_v)
  currUnion_v <- union(currFClones_v, currSClones_v)
  currAbsent_v <- setdiff(totalClones_v, currUnion_v)
  
  ## Calculate Jaccard (intersection / union)
  currJaccard_v <- round(length(currInt_v) / length(currUnion_v), digits = 4)
  
  ## Calculate BUB
  ## Numerator - nIntersect + sqrt(nIntersect * nAbsent)
  ## Denominator - nFirst + nSecond - nIntersect + sqrt(nIntersect * nAbsent) ||||||| aka nUnion + sqrt(nIntersect * nAbsent)
  currNum_BUB <- length(currInt_v) + sqrt(length(currInt_v) * length(currAbsent_v))
  currDenom_BUB <- length(currUnion_v) + sqrt(length(currInt_v) * length(currAbsent_v))
  currBUB_v <- round(currNum_BUB / currDenom_BUB, digits = 4)
  
  ## Get other values
  currFNClones_v <- length(currFClones_v)
  currSNClones_v <- length(currSClones_v)
  currNClones_v <- length(totalClones_v)
  currNInt_v <- length(currInt_v)
  currNUnion_v <- length(currUnion_v)
  
  ## Get row to update
  ## Treat, experiment, mouse will get all of the samples from the pair. These are also defined at the beginning of each
  ## Iteration of totalPairs_dt
  ## Both tissue and sample should be used to differentiate between the samples of the current pair
  ## Can't use just sample because sometimes they overlap. Want to use info from first smaple
  row1_v <- which(output_dt[[animalCol_v]] == currAnimal_v &
                   output_dt[[tissueCol_v]] == currFTissue_v &
                   output_dt[[sampleCol_v]] == currFSample_v)
  row2_v <- which(output_dt[[animalCol_v]] == currAnimal_v &
                    output_dt[[tissueCol_v]] == currSTissue_v &
                    output_dt[[sampleCol_v]] == currSSample_v)
  rows_v <- c(row1_v, row2_v)
        
  ## Update 
  for (k in 1:length(rows_v)){
    row_v <- rows_v[k]
    ## Conditional values
    if (k == 1){
      set(output_dt, i = row_v, j = "PairSamples", value = currSecond_v)
      set(output_dt, i = row_v, j = "ClonesA", value = currFNClones_v)
      set(output_dt, i = row_v, j = "ClonesB", value = currSNClones_v)
    } else {
      set(output_dt, i = row_v, j = "PairSamples", value = currFirst_v)
      set(output_dt, i = row_v, j = "ClonesA", value = currSNClones_v)
      set(output_dt, i = row_v, j = "ClonesB", value = currFNClones_v)
    } # fi
    ## Constant values
    set(output_dt, i = row_v, j = "BUB_Value", value = currBUB_v)
    set(output_dt, i = row_v, j = "J_Value", value = currJaccard_v)
    set(output_dt, i = row_v, j = "TotalClones", value = currNClones_v)
    set(output_dt, i = row_v, j = "Shared", value = currNInt_v)
    set(output_dt, i = row_v, j = "Union", value = currNUnion_v)
  } # for k
  
  ## Venn Diagram
  # vennDir_v <- mkdir(outDir_v, "venn")
  # fileName_v <- file.path(vennDir_v, paste0("mouse-", currAnimal_v, "_S", currFSample_v, "_S", currSSample_v, ".png"))
  # venn_ls <- list()
  # venn_ls[[currFirst_v]] <- currFClones_v
  # venn_ls[[currSecond_v]] <- currSClones_v
  # venn.diagram(venn_ls, file = fileName_v,
  #              main = paste0("Mouse - ", currAnimal_v),
  #              sub = paste0(currFTissue_v, " vs. ", currSTissue_v),
  #              cat.pos = c(0,0))
  
} # for i

# logDir_v <- mkdir(vennDir_v, "logs")
# mvFiles(vennDir_v, logDir_v)

##############
### OUTPUT ###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##############

wkbkOutName_v <- file.path(outDir_v, wkbkName_v)
append_v <- file.exists(wkbkOutName_v)
print("Append?")
print(append_v)

write.xlsx2(x = output_dt, file = path.expand(file.path(outDir_v, wkbkName_v)),
            col.names = T, row.names = F, sheetName = sheetName_v, append = append_v)
