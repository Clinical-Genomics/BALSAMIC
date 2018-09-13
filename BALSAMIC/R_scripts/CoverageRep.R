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
#suppressPackageStartupMessages(library("bit64"))

option_list <- list(
                    make_option(c("-i", "--infile"), type="character",
                                help="Input coverage analysis of Sambamba for exons.", metavar="character"),
                    make_option(c("--genename"), type="character",
                                help="List of gene symbols comma separeted.", metavar="character"),
                    make_option(c("--ensemble"), type="character",
                                help="List of comma separeted gene ensemble ids, not both.", metavar="character"),
                    make_option(c("-t", "--type"), type="character", default="text",
                                help="Output table type format for exon coverage report [default %default].", metavar="character"),
                    make_option(c("--name"), type="character", help="A name for the output table [default %default].",
                                default="Coverage report"),
                    make_option(c("-o", "--outfile"), type="character",
                                help="In case of PDF, output file name [default infile.Coverage.pdf].", metavar="character"),
                    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
                                help="Print some extra output [default %default]")
                    )

'
    %prog [options]
    
    A script to report the exon coverage summary output from Sambamba with for canonical transcripts (longest
transcript) of each gene. In coverage, the following is taken into account: 1. transcript must be protein coding. 2.
Transcript should not have more than one exon with zero coverage. 3. A canonical transcript is a transcript that is
longest and meets criteria 1 and 2.

' -> usage
opt_parser <- OptionParser(usage = usage, option_list=option_list);
arg <- parse_args(opt_parser)
file <- arg$infile
outfile <- arg$outfile

if (is.null(file)){
  print_help(opt_parser)
  stop("An input is required.", call.=FALSE)
}

if (! is.null(arg$genename) & ! is.null(arg$ensemble) ){
  stop("Provide gene or ensemble id, not  both", call.=FALSE)
}

if ( is.null(arg$genename) & is.null(arg$ensemble) ){
  stop("You provide a list of genes or ensemble IDs.", call.=FALSE)
}

if ( arg$verbose ) {
  options( warn = 0 )
} else {
  options( warn = -1 )
}

if ( arg$verbose ) {
  write("Read coverage file.", stderr())
}

sample.coverage = fread(file, showProgress=F)

fragLength = 100
if (! is.null(arg$genename)) {
  genelist = unlist(strsplit(arg$genename, ","))
  dt.gene = sample.coverage[F10 %in% genelist,]
}

if (! is.null(arg$ensemble)) {
  genelist = unlist(strsplit(arg$ensemble, ","))
  dt.gene = sample.coverage[F11 %in% genelist,]
}

dt.gene = dt.gene[,
                  .(exonCount = .N,
                    readPerExon = sum(readCount)/.N,
                    meanExonCoverage = mean(readCount*fragLength/(abs(chromStart-chromEnd))),
                    medianExonCoverage = median(readCount*fragLength/(abs(chromStart-chromEnd))),
                    readPerbpPerExon = sum(readCount*fragLength)/(F7/.N),
                    txID = F3,
                    geneID = F11,
                    txLength = F7,
                    geneName = F10,
                    txType = F6,
                    txStatus = F9,
                    totalRead = sum(readCount),
                    zeroExonCov = sum(readCount==0),
                    zeroExonCovMid = !(any(which(!readCount)==length(readCount)) 
                    || any(which(!readCount)==1)),
                    zeroExonCovLastFirst = any(which(!readCount)==length(readCount)) 
                    || any(which(!readCount)==1)
                   ),
                  by=.(F3, F6, F7, F9, F10, F11)
                 ]

dt.geme = dt.gene[zeroExonCov <= 1 & txType=="protein_coding",]

dt.gene = dt.gene[,
                  .("Gene" = geneName,
                    txID,
                    "tx_exonCount" = paste0(txID, "_", exonCount), 
                    "tx length" = txLength,
                    txLength,
                    maxLength = max(txLength),
                    "tx type" = txType,
                    "tx status" = txStatus,
                    exonCount,
                    "read per exon" = readPerExon,
                    readPerbpPerExon,
                    "Median exon cov" = medianExonCoverage, 
                    meanExonCoverage,
                    totalRead,
                    zeroExonCov,
                    zeroExonCovLastFirst,
                    zeroExonCovMid,
                    "Exon zero cov" = paste0(zeroExonCov,
                                               " (",
                                               as.integer(zeroExonCovLastFirst),
                                               " / ",
                                               as.integer(zeroExonCovMid),
                                               ")"),
                    maxTxReadCount = max(totalRead)
                   ),keyby=geneID]

dt.gene = dt.gene[maxLength==txLength,
                  !c("maxLength",
                     "zeroExonCovLastFirst",
                     "zeroExonCovMid",
                     "txLength",
                     "geneID",
                     "txID",
                     "exonCount",
                     "totalRead",
                     "readPerbpPerExon",
                     "meanExonCoverage",
                     "zeroExonCov",
                     "maxTxReadCount")
                 ]

stargazer(dt.gene, summary = FALSE, type = arg$type, title = arg$name,
          notes = c("Exon zero cov shows the number of exons with zero coverage in total (first or last exon / any intermediate exons)"), 
          table.placement = "H",
          digit.separator = "", rownames = F, style = "io", float = T,
          header = F, out.header = F)
