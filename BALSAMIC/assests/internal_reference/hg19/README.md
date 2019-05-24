Files in this directory are all the reference files used internally by BALSAMIC

### ncbi

```
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ncbiRefSeqPsl.txt.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/seqNcbiRefSeq.txt.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ncbiRefSeq.txt.gz
```

#### seqNcbiRefSeq

header:
`id acc size  gb_date extFile file_offset file_size`

#### ncbiRefSeqPsl

header:
`bin  matches misMatches  repMatches  nCount  qNumInsert  qBaseInsert tNumInsert  tBaseInsert strand  qName qSize qStart
qEnd  tName tSize tStart  tEnd  blockCount  blockSizes  qStarts tStarts`

#### ncbiRefSeq

header:
`bin  name  chrom strand  txStart txEnd cdsStart  cdsEnd  exonCount exonStarts  exonEnds  score name2 cdsStartStat
cdsEndStat  exonFrames`
