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

option_list <- list(
                    make_option(c("-i", "--infile"), type="character",
                                help="Input coverage analysis of Sambamba for exons.", metavar="character"),
                    make_option(c("-o", "--outfile"), type="character",
                                help="In case of PDF, output file name [default infile.Coverage.pdf].", metavar="character"),
                    make_option(c("--avgcov"), type="integer",
                                help="Average coverage of sample. If it's not provided, an average will calculated from input file"),
                    make_option(c("--covline"), type="integer",
                                help="Plot coverage and normalized coverage plot for bed regions in the input file  [default %default]", metavar="character", default=100),
                    make_option(c("--title"), type="character", help="plot title.", metavar="character", default= "Sample"),
                    make_option(c("-f", "--fontsize"), type="integer", default=12,
                                help="Fontsize as an input to pointsize in pdf() for heigh and width [default %default]"),
                    make_option(c("-r", "--resolution"), type="integer", default=7,
                                help="Print image resolution in inches, as an input to pdf() for heigh and width [default %default]"),
                    make_option(c("-v", "--verbose"), action="store_true", default=FALSE,
                                help="Print some extra output [default %default]")
                    )

'
    %prog [options]
    Coverage plot for sambamba depth output file.
' -> usage
opt_parser <- OptionParser(usage = usage, option_list=option_list);
arg <- parse_args(opt_parser)
file <- arg$infile
outfile <- arg$outfile

if (is.null(file)){
  print_help(opt_parser)
  stop("An input is required.", call.=FALSE)
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
covSample = sample.coverage

if (is.null(arg$avgcov) ) {
  arg$avgcov = mean(covSample$meanCoverage)
}

pdf(arg$outfile, width = arg$resolution, height = arg$resolution, pointsize = arg$fontsize)

par(mfrow = c(2,2))

hist(covSample$meanCoverage, breaks = 100, xlab = "Coverage", main = arg$title)
abline(v=arg$covline,col="red")

hist(covSample$meanCoverage/arg$avgcov, breaks = 100, xlab = "Normalized Coverage", main = arg$title)
abline(v=arg$covline/arg$avgcov,col="red")

garbage <- dev.off()
