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
                    make_option(c("--genename"), type="character",
                                help="List of gene symbols comma separeted.", metavar="character"),
                    make_option(c("--ensemble"), type="character",
                                help="List of comma separeted gene ensemble ids, not both.", metavar="character"),
                    make_option(c("-t", "--type"), type="character", default="text",
                                help="Output table type format for exon coverage report [default %default].", metavar="character"),
                    make_option(c("-r", "--resolution"), type="integer", default=7,
                                help="Print image resolution in inches, as an input to pdf() for heigh and width [default %default]"),
                    make_option(c("-f", "--fontsize"), type="integer", default=12,
                                help="Fontsize as an input to pointsize in pdf() for heigh and width [default %default]"),
                    make_option(c("-o", "--outfile"), type="character",
                                help="In case of PDF, output file name [default infile.Coverage.pdf].", metavar="character"),
                    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
                                help="Print some extra output [default %default]")
                    )

'
    %prog [options]
    
    A script to report the exon coverage summary output from Sambamba with for canonical transcripts (longest
transcript) of each gene.

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

sample.coverage = fread(file)


if (! is.null(arg$genename)) {
  genelist = unlist(strsplit(arg$genename, ","))
  dt.gene = sample.coverage[F10 %in% genelist & F9 == "KNOWN" & F6 == "protein_coding", 
                            .(exonCount = .N,
                              readPerbp = sum(readCount)/F7,
                              meanExonCoverage = mean(readCount/(abs(chromStart-chromEnd))),
                              medianExonCoverage = median(readCount/(abs(chromStart-chromEnd))),
                              readPerbpPerExon = sum(readCount)/(F7/.N),
                              totalRead = sum(readCount)),
                            by=.(F3, F7, F11, F10)][, .(txID = F3, txLength = F7, geneName = F10,
                                                        exonCount, readPerbp, readPerbpPerExon, totalRead,
                                                        maxLength = max(F7), maxRead = max(totalRead),
                                                        meanExonCoverage, medianExonCoverage),
                                                    keyby=F11][txLength==maxLength, !c("maxLength")]
}

if (! is.null(arg$ensemble)) {
  genelist = unlist(strsplit(arg$ensemble, ","))
  dt.gene = sample.coverage[F11 %in% genelist & F9 == "KNOWN" & F6 == "protein_coding", 
                            .(exonCount = .N,
                              readPerbp = sum(readCount)/F7,
                              meanExonCoverage = mean(readCount/(abs(chromStart-chromEnd))),
                              medianExonCoverage = median(readCount/(abs(chromStart-chromEnd))),
                              readPerbpPerExon = sum(readCount)/(F7/.N),
                              totalRead = sum(readCount)),
                            by=.(F3, F7, F11, F10)][, .(txID = F3, txLength = F7, geneName = F10,
                                                        exonCount, readPerbp, readPerbpPerExon, totalRead,
                                                        maxLength = max(F7), maxRead = max(totalRead),
                                                        meanExonCoverage, medianExonCoverage),
                                                    keyby=F11][txLength==maxLength, !c("maxLength")]
}
setnames(dt.gene, "F11", "geneID")

stargazer(dt.gene, summary = FALSE, type = "text")

if ( arg$verbose ) {
  write("Converting data frame into mutation", stderr())
}

if ( arg$verbose ) {
  write("Matching signatures with reference.", stderr())
}

#pdf(arg$outfile, width = arg$resolution, height = arg$resolution, pointsize = arg$fontsize)
#plotSignatures(sample_sig)
#garbage <- dev.off()


