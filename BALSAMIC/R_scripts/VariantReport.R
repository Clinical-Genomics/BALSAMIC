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
                    make_option(c("--dp"), type="integer", help="Total read depth filter [default %default].",
                                default=100),
                    make_option(c("-F","--afmax"), type="double", help="Maximum tumor AF filter [default %default].",
                                default=0.05),
                    make_option(c("-f","--afmin"), type="double", help="Minimum tumor AF filter [default %default].",
                                default=0.01),
                    make_option(c("-a","--tumorad"), type="double", help="Allelic depth for alternative allele in tumor
                                [default %default].", default=10),
                    make_option(c("-m","--inMVL"), type="logical", help="Flag to filter variant in MVL or not [default
                                %default].", default=FALSE),
                    make_option(c("--vartype"), type="character", help="Variant type filter. The value must exist in the
                                TYPE column [default %default]", default="SNP"),
                    make_option(c("--varcaller"), type="character", help="Variant caller name filter. Choose from:
                                STRELKA, MUTECT2, VARDICT, or ANY. Use multiple variant caller names sepraterd by comma and no space in between. [default %default].",
                                default="VARDICT"),
                    make_option(c("--ann"), type="character", help="Annotation string to exact match and filter.
                                [default %default].", default="missense_variant"),
                    make_option(c("--name"), type="character", help="A name for the output table [default %default].",
                                default="Variant filter table"),
                    make_option(c("--inExon"), type="logical", help="A flag to select variants that are only found in
                                exons [default %default]", default=TRUE),
                    make_option(c("--inGene"), type="logical", help="A flag to select variants that have a gene symbol
                                annotation [default %default]", default=TRUE),
                    make_option(c("--exportGene"), type="logical", help="A flag to not output the table, instead comma
                                separated list of genes [default %default]", default=FALSE),
                    make_option(c("-t", "--type"), type="character", default="text",
                                help="Output table fortmat type. Choose from: text, latex, html. And output file name is required for html and latex [default %default].", metavar="character"),
                    make_option(c("-o", "--outfile"), type="character",
                                help="In case of PDF, output file name.", metavar="character"),
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

if (is.null(outfile) && arg$type != "text"){
  print_help(opt_parser)
  stop("An output is required for non text output type.", call.=FALSE)
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
       .("Chr:Pos" = paste0(CHROM,":",POS),
         "Ref/Alt" = paste0(REF,"/",ALT),
         "Caller" = ConcatVarCall(paste(unique(c(CALLER)), collapse = "/")),
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
dt = dt[,c("Chr:Pos", "Ref/Alt", "Caller", "DP (Ref/Alt)", "AF", "Gene", "Consequence", "Features")]

if (arg$exportGene) {
  print(paste(unlist(unique(dt[, c("Gene")])), collapse=","))
} else {
  if (arg$type == "text") {

    stargazer(dt, summary = FALSE, type = "text", title = table_name,
              digit.separator = "", rownames = F, style = "io",
              notes = c(paste0("1. A summary of results based on \"", table_name, "\" specification."),
                        paste0("2. Variant callers included: ", tolower(paste(var_caller, collapse = ", ")))))
  } else {
    sink("/dev/null")    
    stargazer(dt, summary = FALSE, type = arg$type, title = table_name,
              digit.separator = "", rownames = F, style = "io", float = T,
              notes = c(paste0("1. A summary of results based on \"", table_name, "\" specification."),
                        paste0("2. Variant callers included: ", tolower(paste(var_caller, collapse = ", ")))),
              header = F, out.header = F, out = arg$outfile)
    sink()
  }
}

