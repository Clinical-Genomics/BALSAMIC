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
suppressPackageStartupMessages(library("deconstructSigs"))

option_list <- list(
                    make_option(c("-i", "--infile"), type="character",
                                help="Input mutation file.", metavar="character"),
                    make_option(c("-o", "--outfile"), type="character",
                                help="output file name.", metavar="character"),
                    make_option(c("-m", "--model"), type="character", default="nature2013",
                                help="Choose between two models: nature2013 or cosmic. Refer to deconstructSigs documentatation. [default %default]"),
                    make_option(c("-r", "--resolution"), type="integer", default=7,
                                help="Print image resolution in inches, as an input to pdf() for heigh and width [default %default]"),
                    make_option(c("-f", "--fontsize"), type="integer", default=12,
                                help="Fontsize as an input to pointsize in pdf() for heigh and width [default %default]"),
                    make_option(c("-s", "--sampleid"), type="integer", default=1,
                                help="Numeric id for sample in input file. [default %default]"),
                    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
                                help="Print extra output [default %default]")
                    )

'
    %prog [options]
    
    A wrapper for deconstructSig to plot mutational signatures.

' -> usage
opt_parser <- OptionParser(usage = usage, option_list=option_list);
arg <- parse_args(opt_parser)
file <- arg$infile
outfile <- arg$outfile 

if (is.null(file) || is.null(outfile) ){
  print_help(opt_parser)
  stop("An input and output files are required.", call.=FALSE)
}

if ( arg$verbose ) {
  options( warn = 0 )
} else {
  options( warn = -1 )
}

if ( arg$verbose ) {
  write("Loading input mutation list into data frame.", stderr())
}

sample.mut.ref = read.table(file , header = T)

if ( arg$verbose ) {
  write("Converting data frame into mutation", stderr())
}

sigs.input <- mut.to.sigs.input(mut.ref = sample.mut.ref,
                                sample.id = "sample",
                                chr = "chrom",
                                pos = "pos",
                                ref = "ref",
                                alt = "alt")

if ( arg$verbose ) {
  write("Matching signatures with reference.", stderr())
}

if ( arg$model == "nature2013" ) {
  sigmodel = signatures.nature2013
} else if ( arg$model == "cosmic" ) {
  sigmodel = signatures.cosmic
} else {
  stop("Unknown model paramters. See help")
}

sample_sig = whichSignatures(tumor.ref = sigs.input,
                             signatures.ref = sigmodel,
                             sample.id = arg$sampleid,
                             contexts.needed = TRUE,
                             tri.counts.method = 'default')

pdf(arg$outfile, width = arg$resolution, height = arg$resolution, pointsize = arg$fontsize)
plotSignatures(sample_sig)
garbage <- dev.off()


