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
          make_option(c("-i", "--infile"), type="character", help="Input coverage analysis of Sambamba for exons.", metavar="character"),
          make_option(c("--mode"), type="character", help="Run mode. Select one from the list: MVL,TMB,VarClass,VarCaller,VarCallerClass [default %default]", default="MVL"),
          make_option(c("--dp"), type="character", help="Total read depth filter [default %default].", default="100"),
          make_option(c("-F","--afmax"), type="character", help="Maximum tumor AF filter [default %default].", default="0.05"),
          make_option(c("-f","--afmin"), type="character", help="Minimum tumor AF filter [default %default].", default="0.01"),
          make_option(c("-a","--tumorad"), type="character", help="Allelic depth for alternative allele in tumor [default %default].", default="10"),
          make_option(c("-m","--inMVL"), type="character", help="Flag to filter variant in MVL or not [default %default].", default="FALSE"),
          make_option(c("--vartype"), type="character", help="Variant type filter. The value must exist in the TYPE column [default %default]", default="SNP"),
          make_option(c("--varcaller"), type="character", help="Variant caller name filter. Choose from: STRELKA, MUTECT2, VARDICT, or ANY. Use multiple variant caller names sepraterd by comma and no space in between. [default %default].", default="VARDICT"),
          make_option(c("--ann"), type="character", help="Annotation string to exact match and filter. [default %default].", default="missense_variant"),
          make_option(c("--name"), type="character", help="A name for the output table [default %default].", default="Variant filter table"),
          make_option(c("--inExon"), type="logical", help="A flag to select variants that are only found in exons [default %default]", default=TRUE),
          make_option(c("--inGene"), type="logical", help="A flag to select variants that have a gene symbol annotation [default %default]", default=TRUE),
          make_option(c("--genomeSize"), type="integer", help="Genome or panel or exome size to calculate TMB"),
          make_option(c("--exclusiveSets"), type="logical", help="A flag to only perform setdiff between different sets of MVL [default %default]", default=FALSE),
          make_option(c("--exportGene"), type="logical", help="A flag to not output the table, instead comma separated list of genes [default %default]", default=FALSE),
          make_option(c("-t", "--type"), type="character", default="text", help="Output table fortmat type. Choose from: text, latex, html. And output file name is required for html and latex [default %default].", metavar="character"),
          make_option(c("-o", "--outfile"), type="character", help="In case of PDF, output file name.", metavar="character"),
          make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="Print some extra output [default %default]")
          )

'
    %prog [options]

    A script to report variants based on series of inputs. Data.table is sometimes faster than Pandas is aggregating
results, thus it was developed in R instead of Python.
' -> usage
opt_parser <- OptionParser(usage = usage, option_list=option_list);
arg <- parse_args(opt_parser)
file <- arg$infile

if (is.null(file)){
  print_help(opt_parser)
  stop("An input is required.", call.=FALSE)
}

if (! arg$verbose ) {
  options( warn = 0 )
} else {
  options( warn = -1 )
}

ConcatVarCall <- function(x) {
    x = gsub("VARDICT","V",x)
    x = gsub("STRELKA","S",x)
    x = gsub("MUTECT2","M",x)
    return(x)
}

trimStr <- function(x) {
    if (nchar(x) > 10) {
        x = paste0(substring(x, 1, 3), "...", substring(x, nchar(x)-3, nchar(x)))
    }
    return(x)
}

sample.coverage = fread(arg$infile, showProgress=F)
sample.coverage[,ID:=paste0(CHROM,"_",POS,"_",REF,"_",ALT)]

dt_excl = data.table()

if (arg$mode == "MVL") {
    var_param = c("afmax","afmin","inMVL","dp","tumorad","name","varcaller","ann","vartype","outfile")

    set_cnt = c()
    for (v in var_param) {
        arg[[v]] = unlist(strsplit(arg[[v]], split=';', fixed=T))
        set_cnt = c(length(unlist(strsplit(arg[[v]],split=';', fixed=T))), set_cnt)
    }

    if (length(unique(set_cnt)) > 1) {
        stop("Number of sets is different among inputs.", call.=FALSE)
    }

    int_vars = c("afmax","afmin","dp","tumorad")
    for (v in int_vars) {
        arg[[v]] = as.numeric(arg[[v]])
    }

    #bool_vars = c("inMVL", "inExon", "inGene")
    bool_vars = c("inMVL")
    for (v in bool_vars) {
        arg[[v]] = as.logical(arg[[v]])
    }
    for (i in 1:unique(set_cnt)) {
        mvl = arg$inMVL[i]
        if (mvl) {
          mvl = "1"
        } else {
          mvl = "."
        }
        dp = arg$dp[i]
        tumor_ad_alt = arg$tumorad[i]
        af_max = arg$afmax[i]
        af_min = arg$afmin[i]
        var_type = unlist(strsplit(arg$vartype[i], ","))
        var_caller = unlist(strsplit(arg$varcaller[i], ","))
        table_name = arg$name[i]
        table_num = arg$num[i]
        outfile = arg$outfile[i]

        var_caller = toupper(var_caller)
        if (any(var_caller %in% "ANY")) {
            var_caller = c("MUTECT2", "VARDICT", "STRELKA")
        }

        dt = sample.coverage[CALLER %in% var_caller
                             & MSK_MVL == mvl
                             & (TUMOR_AD_REF + TUMOR_AD_ALT) >= dp
                             & TUMOR_AD_ALT / (TUMOR_AD_REF + TUMOR_AD_ALT) <= af_max
                             & TUMOR_AD_ALT / (TUMOR_AD_REF + TUMOR_AD_ALT) >= af_min
                             & TUMOR_AD_ALT >= tumor_ad_alt
                             & TYPE %in% var_type]

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

        if (nrow(dt)!=0) {       

            dt = dt[,
                   .("Chr:Pos" = paste0(CHROM,":",POS),
                     "Ref/Alt" = paste0(unlist(lapply(FUN = trimStr, REF)),"/",unlist(lapply(FUN = trimStr, ALT))),
                     "Caller" = ConcatVarCall(paste(unique(c(CALLER)), collapse = "/")),
                     "CallerCount" = length(unique(c(CALLER))),
                     "DP (Ref/Alt)" = paste0(floor(mean(TUMOR_AD_REF + TUMOR_AD_ALT)),
                                             "(",
                                             paste0(floor(mean(TUMOR_AD_REF)),"/", floor(mean(TUMOR_AD_ALT))),
                                             ")"),
                     "AF" = mean(TUMOR_AD_ALT/(TUMOR_AD_REF + TUMOR_AD_ALT)),
                     "Consequence" = paste(unique(c(Consequence)), collapse = ", "),
                     "Protein"=paste(unlist(strsplit(HGVSp,":"))[2], collapse=", "),
                     "Gene" = SYMBOL
                    )
                   ,by=.(ID)]

            dt = unique(dt)
            dt = dt[,c("Chr:Pos", "Ref/Alt", "Caller", "DP (Ref/Alt)", "AF", "Gene", "Consequence", "Protein")]

            if (nrow(dt_excl)==0) {
                dt_excl = dt[0,]
            }
        }


        if (arg$exclusiveSets & unique(set_cnt) > 1 & nrow(dt)>0) {
            dt = fsetdiff(dt, dt_excl, all = FALSE)
            dt_excl = funion(dt, dt_excl)
        }

        if (nrow(dt)==0) {
            write(paste0("No variant were for found for table ", table_name), "")
            write(" ", file = outfile)
        } else {
            if (arg$exportGene) {
                write(paste0("list of genes for table ", table_name, " : ",
                      paste(unlist(unique(dt[, c("Gene")])), collapse=",")), "")
                write(paste(unlist(unique(dt[, c("Gene")])), collapse=","), file = outfile)
            } else {
                if (is.null(arg$outfile) || arg$type == "text") {
                    stargazer(dt, summary = FALSE, type = arg$type, title = table_name,
                              table.placement = "H", digit.separator = "", rownames = F, style = "io", float = T,
                              notes = c(paste0("1. A summary of results based on \"",
                                               table_name, "\" specification."),
                                        paste0("2. Variant callers included: ",
                                               tolower(paste(var_caller, collapse = ", ")))),
                              header = F, out.header = F)
                } else {
                    stargazer(dt, summary = FALSE, title = table_name,
                              table.placement = "H", digit.separator = "", rownames = F, style = "io", float = T,
                              notes = c(paste0("1. A summary of results based on \"",
                                               table_name, "\" specification."),
                                        paste0("2. Variant callers included: ",
                                               tolower(paste(var_caller, collapse = ", ")))),
                              header = F, out.header = F, out = outfile)
                    fwrite(x = dt, file = paste0(outfile, ".csv"))
                }
            }
        }
    }
} else if ( arg$mode == "TMB" ) {
  if (is.null(arg$genomeSize)) {
    stop("Genome/panel size is required.", call.=FALSE)
  }
  genomeLength = arg$genomeSize / 1e6

  var_caller = as.list(unique(sample.coverage[,"CALLER"]))$CALLER
  annotation = c("stop_gained", "stop_lost", "start_lost",
                 "missense_variant", "nonsynonymous_variant",
                 "splice_acceptor_variant", "splice_donor_variant",
                 "splice_donor_5th_base_variant", "splice_site_variant",
                 "splicing_variant", "frameshift_variant")

  dt1 = unique(sample.coverage[CALLER %in% var_caller
                       & Consequence %in% annotation, .(CHROM, POS),
                       by=.(ID, CALLER)])[,.(.N,"TMB"=.N/genomeLength),by=.(CALLER)]

  dt2 = unique(sample.coverage[CALLER %in% var_caller
                       & Consequence %in% annotation, .(CALLER),
                       by=.(ID)][,.(ID)])[,.(.N)]
  dt2[,"CALLER":="ALL"]
  setcolorder(dt2, neworder = c("CALLER", "N"))
  dt = dt2[,.(CALLER, N, "TMB"=N/genomeLength)]

  str_annot = paste(gsub("_", "-", annotation), collapse = ", ")
  dt.TMB = rbind(dt1, dt)
  stargazer(unique(dt.TMB), summary = FALSE, type = arg$type,
            title = "Tumor mutation burden (TMB)",
            digit.separator = "", rownames = F, style = "io",
            header = F, out.header = F, table.placement = "H", float = T, 
            notes = c(paste0("1. Variant callers included: ",
                             tolower(paste(as.list(unique(sample.coverage[,"CALLER"]))$CALLER,
                                           collapse = ", "))),
                      paste0("2. Variant types: ",
                             tolower(paste(as.list(unique(sample.coverage[,"VARIANT_CLASS"]))$VARIANT_CLASS,
                                           collapse=", "))),
                      paste0("3. Only all coding variants (all subchilds of nonsynonymous variants annotation)")))

  if (!is.null(arg$outfile)){
      for (v in var_caller) {
          fwrite(x = unique(sample.coverage[CALLER==v & SYMBOL!=".",
                   .("Chr:Pos" = paste0(CHROM,":",POS),
                     "Ref/Alt" = paste0(unlist(lapply(FUN = trimStr, REF)),"/",unlist(lapply(FUN = trimStr, ALT))),
                     "Caller" = ConcatVarCall(paste(unique(c(CALLER)), collapse = "/")),
                     "CallerCount" = length(unique(c(CALLER))),
                     "DP (Ref/Alt)" = paste0(floor(mean(TUMOR_AD_REF + TUMOR_AD_ALT)),
                                             "(",
                                             paste0(floor(mean(TUMOR_AD_REF)),"/", floor(mean(TUMOR_AD_ALT))),
                                             ")"),
                     "AF" = mean(TUMOR_AD_ALT/(TUMOR_AD_REF + TUMOR_AD_ALT)),
                     "Consequence" = paste(unique(c(Consequence)), collapse = ", "),
                     "Protein" = paste(unlist(unique(HGVSp)), collapse=", "),
                     "Gene" = paste(unlist(unique(SYMBOL)), collapse=", ")
                    )
                   ,by=.(ID)]), file = paste0(arg$outfile, "_", v, ".csv"))
      }
  }
} else if ( arg$mode == "VarClass" ) {
  dt = unique(sample.coverage[,.(ID,CALLER,VARIANT_CLASS)])
  dt.typevars = dt[,.("CALLERCOUNT"=length(unique(c(CALLER))),.N),
                   by=.("CALLERS"=CALLER,VARIANT_CLASS)][order(CALLERS,-VARIANT_CLASS)]
  
  dt = unique(sample.coverage[,.(ID,CALLER,VARIANT_CLASS)])
  dt = dt[,.("CALLERS"="ALL","CALLERCOUNT"=length(unique(c(CALLER))),.N),
          by=(VARIANT_CLASS)]
  
  setcolorder(dt, neworder = c("CALLERS", "VARIANT_CLASS", "CALLERCOUNT","N"))
  dt.typecensus = rbind(dt, dt.typevars)
  
  dt = unique(sample.coverage[,.(ID,CALLER)])
  dt.allcallers = dt[,.("VARIANT_CLASS"="All_types","CALLERCOUNT"=1,.N), by=.("CALLERS"=CALLER) ][order(CALLERS)]
  
  dt.callercensus = rbind(dt.allcallers,dt.typecensus)[order(CALLERS,VARIANT_CLASS,CALLERCOUNT)]
  stargazer(unique(dt.callercensus), summary = FALSE, type = arg$type,
            title = "Variant class summary",
            digit.separator = "", rownames = F, style = "io",
            header = F, out.header = F, table.placement = "H", float = T, 
            notes = c(paste0("1. A summary of variant classes devided by variant class and variant caller"),
                      paste0("2. Variant callers included: ",
                             tolower(paste(as.list(unique(sample.coverage[,"CALLER"]))$CALLER,
                                           collapse = ", "))),
                      paste0("3. Variant types: ",
                             tolower(paste(as.list(unique(sample.coverage[,"VARIANT_CLASS"]))$VARIANT_CLASS,
                                           collapse=", ")))))
} else if ( arg$mode == "VarCaller" ) {
  
  var_caller = as.list(unique(sample.coverage[,"CALLER"]))$CALLER
  var_caller_combn = do.call("c", lapply(seq_along(var_caller),
                                         function(i) {combn(var_caller, i, simplify = F)}))

  dt = unique(sample.coverage[,.(ID,CALLER,VARIANT_CLASS)])
  dt = dt[,.("CALLERS"=paste(unique(c(CALLER)),collapse = "-"),
                          "CALLERCOUNT"=length(unique(c(CALLER)))),
                       by=.(ID)][,.(.N),
                                 by=.(CALLERS,CALLERCOUNT)
                                ]

  dt.venn = dt[order(CALLERCOUNT,CALLERS)]
  stargazer(unique(dt.venn), summary = FALSE, type = arg$type,
            title = "Variant caller summary",
            digit.separator = "", rownames = F, style = "io",
            header = F, out.header = F, table.placement = "H", float = T, 
            notes = c(paste0("1. A summary of exclusive variant types ",
                             "devided by variant callers that identified them"),
                      paste0("2. Variant callers included: ",
                             tolower(paste(as.list(unique(sample.coverage[,"CALLER"]))$CALLER,
                                           collapse = ", "))),
                      paste0("3. Variant types: ",
                             tolower(paste(as.list(unique(sample.coverage[,"VARIANT_CLASS"]))$VARIANT_CLASS,
                                           collapse=", ")))))
} else if ( arg$mode == "VarCallerClass" ){
  
  var_caller = as.list(unique(sample.coverage[,"CALLER"]))$CALLER
  var_caller_combn = do.call("c", lapply(seq_along(var_caller),
                                         function(i) {combn(var_caller, i, simplify = F)}))

  dt = unique(sample.coverage[,.(ID,CALLER,VARIANT_CLASS)])
  dt = dt[,.("CALLERS"=paste(unique(c(CALLER)),collapse = "-"),
                          "CALLERCOUNT"=length(unique(c(CALLER))),
                          VARIANT_CLASS),
                       by=.(ID)][,.(.N),
                                 by=.(CALLERS,VARIANT_CLASS,CALLERCOUNT)
                                ]

  dt.venn = dt[order(CALLERCOUNT,CALLERS,VARIANT_CLASS)]
  stargazer(unique(dt.venn), summary = FALSE, type = arg$type,
            title = "Variant caller summary by class",
            digit.separator = "", rownames = F, style = "io",
            header = F, out.header = F, table.placement = "H", float = T, 
            notes = c(paste0("1. A summary of exclusive variant types ",
                             "devided by variant callers that identified them"),
                      paste0("2. Variant callers included: ",
                             tolower(paste(as.list(unique(sample.coverage[,"CALLER"]))$CALLER,
                                           collapse = ", "))),
                      paste0("3. Variant types: ",
                             tolower(paste(as.list(unique(sample.coverage[,"VARIANT_CLASS"]))$VARIANT_CLASS,
                                           collapse=", ")))))
} else {
  stop("Run mode not recognized", call.=FALSE) 
}
