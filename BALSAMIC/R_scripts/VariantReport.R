#!/usr/bin/env Rscript
#  Copyright 2018, Hassan Foroughi Asl <hassan.foroughi@scilifelab.se>
#  
#  This file is free software: you may copy, redistribute and/or modify it  
#  under the terms of the GNU General Public License as published by the  
#  Free Software Foundation, either version 2 of the License, or (at your  
#  option) any later version.  
#  
#  This file is distributed in the hope that it will be useful, but  
#  WITHOUT ANY WARRANTY; without even the implied warranty of  
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU  
#  General Public License for more details.  
#  
#  You should have received a copy of the GNU General Public License  
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.  

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("stargazer"))

option_list <- list(
                    make_option(c("-i", "--infile"), type="character",
                                help="Input coverage analysis of Sambamba for exons.", metavar="character"),
                    make_option(c("--dp"), type="integer"),
                    make_option(c("-F","--afmax"), type="double"),
                    make_option(c("-f","--afmin"), type="double"),
                    make_option(c("-a","--tumorad"), type="double"),
                    make_option(c("-m","--inMVL"), type="logical"),
                    make_option(c("--vartype"), type="character"),
                    make_option(c("--varcaller"), type="character"),
                    make_option(c("--ann"), type="character"),
                    make_option(c("--name"), type="character"),
                    make_option(c("--num"), type="integer"),
                    make_option(c("--inExon"), type="logical"),
                    make_option(c("--inGene"), type="logical"),
                    make_option(c("--exportGene"), type="logical"),
                    make_option(c("-t", "--type"), type="character", default="text",
                                help="Output table type format for exon coverage report [default %default].", metavar="character"),
                    make_option(c("-o", "--outfile"), type="character",
                                help="In case of PDF, output file name [default infile.Coverage.pdf].", metavar="character"),
                    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
                                help="Print some extra output [default %default]")
                    )

'
    %prog [options]
    
    A script to report variants based on series of inputs. Data.table is sometimes faster than Pandas is aggregating
results, thus it was developed in R instead of Python.

' -> usage
opt_parser <- OptionParser(usage = usage, option_list=option_list);
arg <- parse_args(opt_parser)
file <- arg$infile
outfile <- arg$outfile

if (is.null(file)){
  print_help(opt_parser)
  stop("An input is required.", call.=FALSE)
}

if (! arg$verbose ) {
  options( warn = 0 )
} else {
  options( warn = -1 )
}

dp = arg$dp
tumor_ad_alt = arg$tumorad
if (arg$inMVL) {
  mvl = "1"
} else {
  mvl = "."
}
af_max = arg$afmax
af_min = arg$afmin
var_type = arg$vartype
#c("stop_gained", "stop_lost", "start_lost",
#                  "missense_variant", "nonsynonymous_variant",
#                  "splice_acceptor_variant", "splice_donor_variant",
#                  "splice_donor_5th_base_variant", "splice_site_variant", "splicing_variant")
#annotation_indel = c("frameshift_variant", "frameshift", "non-frameshift")
var_caller = unlist(strsplit(arg$varcaller, ","))
table_name = arg$name
table_num = arg$num

var_caller = toupper(var_caller)
if (any(var_caller %in% "ANY")) {
    var_caller = c("MUTECT2", "VARDICT", "STRELKA")
}

ConcatVarCall <- function(x) {
    x = gsub("VARDICT","V",x)
    x = gsub("STRELKA","S",x)
    x = gsub("MUTECT2","M",x)
    return(x)
}

sample.coverage = fread(arg$infile)
sample.coverage[,ID:=paste0(CHROM,"_",POS,"_",REF,"_",ALT)]
options(repr.matrix.max.cols = 80)

dt = sample.coverage[CALLER %in% var_caller
                     & MSK_MVL == mvl
                     & (TUMOR_AD_REF + TUMOR_AD_ALT) >= dp
                     & TUMOR_AD_ALT / (TUMOR_AD_REF + TUMOR_AD_ALT) <= af_max
                     & TUMOR_AD_ALT / (TUMOR_AD_REF + TUMOR_AD_ALT) >= af_min
                     & TUMOR_AD_ALT >= tumor_ad_alt
                     & TYPE == var_type]

if (! is.null(arg$ann)) {
  var_ann = unlist(strsplit(arg$ann, ",")) 
  dt = dt[Consequence %in% var_ann]
}

if (arg$inExon) {
  dt = dt[EXON != "."]
}

if (arg$inGene) {
  dt = dt[SYMBOL != "."]
}

dt = dt[,
       .(Chr = CHROM,
         Pos = POS,
         "Ref/Alt" = paste0(REF,"/",ALT),
         "Var Caller" = ConcatVarCall(paste(unique(c(CALLER)), collapse = "/")),
         "CallerCount" = length(unique(c(CALLER))),
         "DP (Ref/Alt)" = paste0(floor(mean(TUMOR_AD_REF + TUMOR_AD_ALT)),
                                 "(",
                                 paste0(floor(mean(TUMOR_AD_REF)),"/", floor(mean(TUMOR_AD_ALT))),
                                 ")"),
         "AF" = mean(TUMOR_AD_ALT/(TUMOR_AD_REF + TUMOR_AD_ALT)),
         "Consequence" = paste(unique(c(Consequence)), collapse = ", "),
         "Gene" = SYMBOL,
         "Features" = paste(unique(c(BIOTYPE)), collapse = ", ")
        )
       ,by=.(ID)]

dt = unique(dt)
dt = dt[,c("Chr", "Pos", "Ref/Alt", "Var Caller", "DP (Ref/Alt)", "AF", "Gene", "Consequence", "Features")]

if (arg$exportGene) {
  print(paste(unlist(unique(dt[, c("Gene")])), collapse=","))
}
stargazer(dt, summary = FALSE, type = "text", title = paste0("Table ", table_num, ": ", table_name),
          digit.separator = "",
          notes = c(paste0("1. A summary of results based on \"", table_name, "\" specification."),
                    paste0("2. Variant callers included: ", tolower(paste(var_caller, collapse = ", ")))))#,
#          out.header = T, out = "notebook_data/OTT3801A3_R.sorted.rmdup.exon.cov.bed.report")

