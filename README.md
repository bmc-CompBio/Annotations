
## Bioconductor Resources

https://www.bioconductor.org/help/workflows/annotation-data/

https://bioconductor.org/help/workflows/annotation/Annotation_Resources/

## Libraries


```r
library(TxDb.Dmelanogaster.UCSC.dm6.ensGene)
library(org.Dm.eg.db)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(GO.db)
library(rtracklayer)
```

## TxDb


```r
# TxDb object

txdb <- TxDb.Dmelanogaster.UCSC.dm6.ensGene
```


```
## TxDb object:
## # Db type: TxDb
## # Supporting package: GenomicFeatures
## # Data source: UCSC
## # Genome: dm6
## # Organism: Drosophila melanogaster
## # Taxonomy ID: 7227
## # UCSC Table: ensGene
## # UCSC Track: Ensembl Genes
## # Resource URL: http://genome.ucsc.edu/
## # Type of Gene ID: Ensembl gene ID
## # Full dataset: yes
## # miRBase build ID: NA
## # transcript_nrow: 34703
## # exon_nrow: 83433
## # cds_nrow: 62685
## # Db created by: GenomicFeatures package from Bioconductor
## # Creation time: 2017-04-04 22:11:07 +0000 (Tue, 04 Apr 2017)
## # GenomicFeatures version at creation time: 1.27.13
## # RSQLite version at creation time: 1.1-2
## # DBSCHEMAVERSION: 1.1
```



```r
# extract ranges

fly_genes <- genes(txdb)

fly_exons <- exonsBy(txdb, by = "gene")

fly_fiveUTR <- fiveUTRsByTranscript(txdb)
```



```
## GRanges object with 17522 ranges and 1 metadata column:
##               seqnames               ranges strand |     gene_id
##                  <Rle>            <IRanges>  <Rle> | <character>
##   FBgn0000003    chr3R [ 6822498,  6822796]      + | FBgn0000003
##   FBgn0000008    chr2R [22136968, 22172834]      + | FBgn0000008
##   FBgn0000014    chr3R [16807214, 16830049]      - | FBgn0000014
##   FBgn0000015    chr3R [16927210, 16972236]      - | FBgn0000015
##   FBgn0000017    chr3L [16615866, 16647882]      - | FBgn0000017
##           ...      ...                  ...    ... .         ...
##   FBgn0267791    chr2R [20587551, 20606854]      + | FBgn0267791
##   FBgn0267792    chr2R [ 8600415,  8612295]      + | FBgn0267792
##   FBgn0267793    chr2R [16297449, 16306695]      - | FBgn0267793
##   FBgn0267794    chr3L [17934032, 17946114]      + | FBgn0267794
##   FBgn0267795    chr3L [14053215, 14071583]      + | FBgn0267795
##   -------
##   seqinfo: 1870 sequences (1 circular) from dm6 genome
```

```
## GRangesList object of length 17522:
## $FBgn0000003
## GRanges object with 1 range and 2 metadata columns:
##       seqnames             ranges strand |   exon_id   exon_name
##          <Rle>          <IRanges>  <Rle> | <integer> <character>
##   [1]    chr3R [6822498, 6822796]      + |     49264        <NA>
##
## $FBgn0000008
## GRanges object with 13 ranges and 2 metadata columns:
##        seqnames               ranges strand | exon_id exon_name
##    [1]    chr2R [22136968, 22137026]      + |   22500      <NA>
##    [2]    chr2R [22136968, 22137208]      + |   22501      <NA>
##    [3]    chr2R [22137433, 22138251]      + |   22502      <NA>
##    [4]    chr2R [22138000, 22138251]      + |   22503      <NA>
##    [5]    chr2R [22151654, 22151695]      + |   22508      <NA>
##    ...      ...                  ...    ... .     ...       ...
##    [9]    chr2R [22170778, 22171259]      + |   22512      <NA>
##   [10]    chr2R [22170778, 22171985]      + |   22513      <NA>
##   [11]    chr2R [22172082, 22172252]      + |   22514      <NA>
##   [12]    chr2R [22172316, 22172433]      + |   22515      <NA>
##   [13]    chr2R [22172497, 22172834]      + |   22516      <NA>
##
## ...
## <17520 more elements>
## -------
## seqinfo: 1870 sequences (1 circular) from dm6 genome
```

```
## GRangesList object of length 29726:
## $1
## GRanges object with 1 range and 3 metadata columns:
##       seqnames       ranges strand |   exon_id   exon_name exon_rank
##          <Rle>    <IRanges>  <Rle> | <integer> <character> <integer>
##   [1]    chr2L [7529, 7679]      + |         1        <NA>         1
##
## $2
## GRanges object with 1 range and 3 metadata columns:
##       seqnames       ranges strand | exon_id exon_name exon_rank
##   [1]    chr2L [7529, 7679]      + |       1      <NA>         1
##
## $3
## GRanges object with 1 range and 3 metadata columns:
##       seqnames       ranges strand | exon_id exon_name exon_rank
##   [1]    chr2L [7529, 7679]      + |       1      <NA>         1
##
## ...
## <29723 more elements>
## -------
## seqinfo: 1870 sequences (1 circular) from dm6 genome
```


## OrgDb



```r
# OrgDb object

orgdb <- org.Dm.eg.db
```


```
## OrgDb object:
## | DBSCHEMAVERSION: 2.1
## | Db type: OrgDb
## | Supporting package: AnnotationDbi
## | DBSCHEMA: FLY_DB
## | ORGANISM: Drosophila melanogaster
## | SPECIES: Fly
## | EGSOURCEDATE: 2017-Nov6
## | EGSOURCENAME: Entrez Gene
## | EGSOURCEURL: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA
## | CENTRALID: EG
## | TAXID: 7227
## | GOSOURCENAME: Gene Ontology
## | GOSOURCEURL: ftp://ftp.geneontology.org/pub/go/godatabase/archive/latest-lite/
## | GOSOURCEDATE: 2017-Nov01
## | GOEGSOURCEDATE: 2017-Nov6
## | GOEGSOURCENAME: Entrez Gene
## | GOEGSOURCEURL: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA
## | KEGGSOURCENAME: KEGG GENOME
## | KEGGSOURCEURL: ftp://ftp.genome.jp/pub/kegg/genomes
## | KEGGSOURCEDATE: 2011-Mar15
## | GPSOURCENAME: UCSC Genome Bioinformatics (Drosophila melanogaster)
## | GPSOURCEURL:
## | GPSOURCEDATE: 2014-Dec12
## | FBSOURCEDATE: 2017-Feb01
## | FBSOURCENAME: Flybase
## | FBSOURCEURL: ftp://ftp.flybase.net/releases/current/precomputed_files/genes/
## | ENSOURCEDATE: 2017-Aug23
## | ENSOURCENAME: Ensembl
## | ENSOURCEURL: ftp://ftp.ensembl.org/pub/current_fasta
```


```r
# columns

columns(orgdb)
```

```
##  [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "FLYBASE"      "FLYBASECG"    "FLYBASEPROT"  "GENENAME"     "GO"           "GOALL"        "MAP"          "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PMID"         "REFSEQ"       "SYMBOL"       "UNIGENE"      "UNIPROT"
```


```r
# mapping gene ids

favorite_gene <- "roX2"

mapIds(orgdb,
       keys = favorite_gene,
       keytype = "SYMBOL",
       column = "FLYBASE",
       multiVals = "first")
```

```
##          roX2
## "FBgn0019660"
```


```r
# update gene ranges

fly_genes$symbol <-  mapIds(orgdb,
                            keys = fly_genes$gene_id,
                            keytype = "FLYBASE",
                            column = "SYMBOL",
                            multiVals = "first")
```



```
## GRanges object with 17522 ranges and 2 metadata columns:
##               seqnames               ranges strand |     gene_id         symbol
##                  <Rle>            <IRanges>  <Rle> | <character>    <character>
##   FBgn0000003    chr3R [ 6822498,  6822796]      + | FBgn0000003 7SLRNA:CR32864
##   FBgn0000008    chr2R [22136968, 22172834]      + | FBgn0000008              a
##   FBgn0000014    chr3R [16807214, 16830049]      - | FBgn0000014          abd-A
##   FBgn0000015    chr3R [16927210, 16972236]      - | FBgn0000015          Abd-B
##   FBgn0000017    chr3L [16615866, 16647882]      - | FBgn0000017            Abl
##           ...      ...                  ...    ... .         ...            ...
##   FBgn0267791    chr2R [20587551, 20606854]      + | FBgn0267791        HnRNP-K
##   FBgn0267792    chr2R [ 8600415,  8612295]      + | FBgn0267792            rgr
##   FBgn0267793    chr2R [16297449, 16306695]      - | FBgn0267793        CR45232
##   FBgn0267794    chr3L [17934032, 17946114]      + | FBgn0267794        CR43174
##   FBgn0267795    chr3L [14053215, 14071583]      + | FBgn0267795            Frl
##   -------
##   seqinfo: 1870 sequences (1 circular) from dm6 genome
```


```r
# GO

favorite_gene_GO <- mapIds(orgdb,
                           keys = favorite_gene,
                           keytype = "SYMBOL",
                           column = "GO",
                           multiVals = "list")
```


```
## $roX2
##  [1] "GO:0000805" "GO:0000805" "GO:0005515" "GO:0007549" "GO:0007549" "GO:0009047" "GO:0016456" "GO:0016457" "GO:0031936" "GO:0035613" "GO:0072487"
```

## GODb



```r
# GODb object

GO.db
```

```
## GODb object:
## | GOSOURCENAME: Gene Ontology
## | GOSOURCEURL: ftp://ftp.geneontology.org/pub/go/godatabase/archive/latest-lite/
## | GOSOURCEDATE: 2017-Nov01
## | Db type: GODb
## | package: AnnotationDbi
## | DBSCHEMA: GO_DB
## | GOEGSOURCEDATE: 2017-Nov6
## | GOEGSOURCENAME: Entrez Gene
## | GOEGSOURCEURL: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA
## | DBSCHEMAVERSION: 2.1
```


```r
# convert ids

select(GO.db, favorite_gene_GO[[1]], "TERM", "GOID")
```

```
##          GOID                                                                                                    TERM
## 1  GO:0000805                                                                                            X chromosome
## 2  GO:0000805                                                                                            X chromosome
## 3  GO:0005515                                                                                         protein binding
## 4  GO:0007549                                                                                     dosage compensation
## 5  GO:0007549                                                                                     dosage compensation
## 6  GO:0009047                                                  dosage compensation by hyperactivation of X chromosome
## 7  GO:0016456                              X chromosome located dosage compensation complex, transcription activating
## 8  GO:0016457 dosage compensation complex assembly involved in dosage compensation by hyperactivation of X chromosome
## 9  GO:0031936                                                              negative regulation of chromatin silencing
## 10 GO:0035613                                                                                   RNA stem-loop binding
## 11 GO:0072487                                                                                             MSL complex
```

## BSgenome



```r
# genome objects

fly_genome <- BSgenome.Dmelanogaster.UCSC.dm6
```



```
##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
##  Fly genome:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
##  # organism: Drosophila melanogaster (Fly)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
##  # provider: UCSC                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
##  # provider version: dm6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
##  # release date: Aug. 2014                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
##  # release name: BDGP Release 6 + ISO1 MT                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
##  # 1870 sequences:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
##  #   chr2L                  chr2R                  chr3L                  chr3R                  chr4                   chrX                   chrY                   chrM                   chrX_CP007103v1_random chrX_CP007104v1_random chrX_DS483648v1_random chrX_DS483655v1_random chrX_DS483660v1_random chrX_DS483665v1_random chrX_DS483666v1_random chrX_DS483669v1_random chrX_DS483685v1_random chrX_DS483698v1_random chrX_DS483745v1_random chrX_DS483784v1_random chrX_DS483789v1_random chrX_DS483795v1_random chrX_DS483803v1_random chrX_DS483809v1_random chrX_DS483818v1_random chrX_DS483821v1_random chrX_DS483843v1_random chrX_DS483851v1_random chrX_DS483885v1_random chrX_DS483888v1_random chrX_DS483892v1_random chrX_DS483893v1_random chrX_DS483897v1_random chrX_DS483903v1_random chrX_DS483905v1_random chrX_DS483907v1_random chrX_DS483909v1_random chrX_DS483923v1_random chrX_DS483926v1_random chrX_DS483928v1_random chrX_DS483946v1_random chrX_DS483948v1_random chrX_DS483950v1_random
```



```r
# seqinfo

seqinfo(fly_genome)
```

```
## Seqinfo object with 1870 sequences (1 circular) from dm6 genome:
##   seqnames         seqlengths isCircular genome
##   chr2L              23513712      FALSE    dm6
##   chr2R              25286936      FALSE    dm6
##   chr3L              28110227      FALSE    dm6
##   chr3R              32079331      FALSE    dm6
##   chr4                1348131      FALSE    dm6
##   ...                     ...        ...    ...
##   chrUn_DS485998v1       1003      FALSE    dm6
##   chrUn_DS486002v1       1001      FALSE    dm6
##   chrUn_DS486004v1       1001      FALSE    dm6
##   chrUn_DS486005v1       1001      FALSE    dm6
##   chrUn_DS486008v1       1001      FALSE    dm6
```



```r
# extract sequemce

favorite_gene_ranges <- fly_genes[grep(favorite_gene, fly_genes$symbol)]

getSeq(fly_genome, favorite_gene_ranges)
```

```
##   A DNAStringSet instance of length 1
##     width seq                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       names               
## [1]  1368 ATTCGCGGCCTGGTCACACTAAGCTAGGGCTACTTTTTATATCATAAGTCGAGCGTTTAGGTAAGCGAAACAAAAAGAGCTTTAGTTAGAGGTATTAGTTTGGTTGCTATTATTTCTAAATTGAAAACTTCACTTCCTTTGATCCAAAAGACACTGAAAAGACACGTTTGCAGTTGAGTTCATTATTTTCTGGTATACATACACACTGATTACTCATTCAATTGGCATTTTGCTCTTGTTTTTCTCCGATTGCCTTGCACTCGCATAAATTTAACACAAAAAAGAAGTTCGGGGTGTTTAGAATCCATCCACTTGGTACAGTTCCCATCGAGCTGGTGAGTACTCCGCGCAGTGCAACGTATACACACTTGAACTCGAATTCTGGCAGCAAAATGGTGGTGACCTATTAATAACACCCGCTCAATTTTCCCTCGCTTTTCCATTCAGTCTCCAACTAGCTGCAAGTGAGTACGTTACTTCCGA...ACAATACAATACAAGACAAAAAAATGTGTCTTGGAACGCAACATTGTACAAGTCGCAATGCAAACTGAAGTCTTAAAAGACGTGTAAAATGTTGCAAATTAAGCAAATATATATGCATATATGGGTAACGTTTTACGCGCCTTAACCAGTCAAAATACAAAATAAATTGGTAAATTTCATATAACTAGTGAAATGTTATACGAAACTTAACAATTGCCAAATAATACAGATCGATTTAGAGCGAGATGACAATAGAGAGGCGATCTCTCTCGTATACGAGTCTTGAAAAGAAAGAGAAGGCGAACGGTGCTGGCTTAGAGAGAGATGGCAATACTAATTAACTGCAAATACATTTCCGCCATTTTGTTGGCGCTAAAAGTAACGGAAATTCGAGATGCTTTTAGGGCTGCCACCTTGGTTTCCAGGGTGACCAGAACTGACACTTAAATTACATATGACAAATAAAGACTTATCTGCTATTGC FBgn0019660
```


## External Annotations



```r
# TxDb from GTF

txdb <- makeTxDbFromGFF("dmel-all-r6.17.gtf", format = "gtf")
```


```
## TxDb object:
## # Db type: TxDb
## # Supporting package: GenomicFeatures
## # Data source: dmel-all-r6.17.gtf
## # Organism: NA
## # Taxonomy ID: NA
## # miRBase build ID: NA
## # Genome: NA
## # transcript_nrow: 34990
## # exon_nrow: 83728
## # cds_nrow: 62778
## # Db created by: GenomicFeatures package from Bioconductor
## # Creation time: 2018-02-23 13:19:27 +0100 (Fri, 23 Feb 2018)
## # GenomicFeatures version at creation time: 1.30.3
## # RSQLite version at creation time: 2.0
## # DBSCHEMAVERSION: 1.2
```


```r
# extract genes

genes(txdb)
```

```
## GRanges object with 17736 ranges and 1 metadata column:
##               seqnames               ranges strand |     gene_id
##                  <Rle>            <IRanges>  <Rle> | <character>
##   FBgn0000003       3R [ 6822498,  6822796]      + | FBgn0000003
##   FBgn0000008       2R [22136968, 22172834]      + | FBgn0000008
##   FBgn0000014       3R [16807214, 16830049]      - | FBgn0000014
##   FBgn0000015       3R [16927212, 16972236]      - | FBgn0000015
##   FBgn0000017       3L [16615866, 16647882]      - | FBgn0000017
##           ...      ...                  ...    ... .         ...
##   FBgn0285994       3L [21576460, 21576661]      - | FBgn0285994
##   FBgn0286004        X [19558607, 19558693]      + | FBgn0286004
##   FBgn0286005       2R [18180151, 18180205]      - | FBgn0286005
##   FBgn0286006       2R [12170383, 12170449]      - | FBgn0286006
##   FBgn0286007       2R [12170163, 12170230]      - | FBgn0286007
##   -------
##   seqinfo: 25 sequences from an unspecified genome; no seqlengths
```


```r
# GRanges from GTF

fly_gtf <- import.gff("dmel-all-r6.17.gtf")
```



```
## GRanges object with 541808 ranges and 9 metadata columns:
##            seqnames               ranges strand |   source        type     score     phase     gene_id gene_symbol transcript_id transcript_symbol           #
##               <Rle>            <IRanges>  <Rle> | <factor>    <factor> <numeric> <integer> <character> <character>   <character>       <character> <character>
##        [1]        X [19961297, 19969323]      + |  FlyBase        gene      <NA>      <NA> FBgn0031081        Nep3          <NA>              <NA>        <NA>
##        [2]        X [19961689, 19968479]      + |  FlyBase        mRNA      <NA>      <NA> FBgn0031081        Nep3   FBtr0070000           Nep3-RA        <NA>
##        [3]        X [19961689, 19961845]      + |  FlyBase        5UTR      <NA>      <NA> FBgn0031081        Nep3   FBtr0070000           Nep3-RA        <NA>
##        [4]        X [19961689, 19961845]      + |  FlyBase        exon      <NA>      <NA> FBgn0031081        Nep3   FBtr0070000           Nep3-RA        <NA>
##        [5]        X [19963955, 19964071]      + |  FlyBase        exon      <NA>      <NA> FBgn0031081        Nep3   FBtr0070000           Nep3-RA        <NA>
##        ...      ...                  ...    ... .      ...         ...       ...       ...         ...         ...           ...               ...         ...
##   [541804]       2L     [779173, 779655]      + |  FlyBase        exon      <NA>      <NA> FBgn0031277     CG13947   FBtr0078005        CG13947-RA        <NA>
##   [541805]       2L     [779204, 779206]      + |  FlyBase start_codon      <NA>         0 FBgn0031277     CG13947   FBtr0078005        CG13947-RA        <NA>
##   [541806]       2L     [779204, 779560]      + |  FlyBase         CDS      <NA>         0 FBgn0031277     CG13947   FBtr0078005        CG13947-RA        <NA>
##   [541807]       2L     [779561, 779563]      + |  FlyBase  stop_codon      <NA>         0 FBgn0031277     CG13947   FBtr0078005        CG13947-RA        <NA>
##   [541808]       2L     [779564, 779655]      + |  FlyBase        3UTR      <NA>      <NA> FBgn0031277     CG13947   FBtr0078005        CG13947-RA        <NA>
##   -------
##   seqinfo: 25 sequences from an unspecified genome; no seqlengths
```



```r
# extract by type

fly_genes <- fly_gtf[fly_gtf$type == "gene"]

fly_cds <- fly_gtf[fly_gtf$type == "CDS"]
```


```
## GRanges object with 17737 ranges and 9 metadata columns:
##           seqnames               ranges strand |   source     type     score     phase     gene_id      gene_symbol transcript_id transcript_symbol           #
##              <Rle>            <IRanges>  <Rle> | <factor> <factor> <numeric> <integer> <character>      <character>   <character>       <character> <character>
##       [1]        X [19961297, 19969323]      + |  FlyBase     gene      <NA>      <NA> FBgn0031081             Nep3          <NA>              <NA>        <NA>
##       [2]        X [20025099, 20025170]      + |  FlyBase     gene      <NA>      <NA> FBgn0052826 tRNA:Pro-CGG-1-1          <NA>              <NA>        <NA>
##       [3]        X [20051294, 20052519]      + |  FlyBase     gene      <NA>      <NA> FBgn0031085           CG9570          <NA>              <NA>        <NA>
##       [4]        X [20094398, 20095767]      + |  FlyBase     gene      <NA>      <NA> FBgn0062565            Or19b          <NA>              <NA>        <NA>
##       [5]        X [20133579, 20138878]      + |  FlyBase     gene      <NA>      <NA> FBgn0031088          CG15322          <NA>              <NA>        <NA>
##       ...      ...                  ...    ... .      ...      ...       ...       ...         ...              ...           ...               ...         ...
##   [17733]       2L     [728362, 730564]      + |  FlyBase     gene      <NA>      <NA> FBgn0011244           Hsp60B          <NA>              <NA>        <NA>
##   [17734]       2L     [749945, 762399]      + |  FlyBase     gene      <NA>      <NA> FBgn0031275        GABA-B-R3          <NA>              <NA>        <NA>
##   [17735]       2L     [773547, 774017]      + |  FlyBase     gene      <NA>      <NA> FBgn0031276          CG12506          <NA>              <NA>        <NA>
##   [17736]       2L     [776536, 776941]      + |  FlyBase     gene      <NA>      <NA> FBgn0040725          CG13946          <NA>              <NA>        <NA>
##   [17737]       2L     [779173, 779655]      + |  FlyBase     gene      <NA>      <NA> FBgn0031277          CG13947          <NA>              <NA>        <NA>
##   -------
##   seqinfo: 25 sequences from an unspecified genome; no seqlengths
```

```
## GRanges object with 160894 ranges and 9 metadata columns:
##            seqnames               ranges strand |   source     type     score     phase     gene_id gene_symbol transcript_id transcript_symbol           #
##               <Rle>            <IRanges>  <Rle> | <factor> <factor> <numeric> <integer> <character> <character>   <character>       <character> <character>
##        [1]        X [19963955, 19964071]      + |  FlyBase      CDS      <NA>         0 FBgn0031081        Nep3   FBtr0070000           Nep3-RA        <NA>
##        [2]        X [19964782, 19964944]      + |  FlyBase      CDS      <NA>         0 FBgn0031081        Nep3   FBtr0070000           Nep3-RA        <NA>
##        [3]        X [19965006, 19965126]      + |  FlyBase      CDS      <NA>         2 FBgn0031081        Nep3   FBtr0070000           Nep3-RA        <NA>
##        [4]        X [19965197, 19965511]      + |  FlyBase      CDS      <NA>         1 FBgn0031081        Nep3   FBtr0070000           Nep3-RA        <NA>
##        [5]        X [19965577, 19966071]      + |  FlyBase      CDS      <NA>         1 FBgn0031081        Nep3   FBtr0070000           Nep3-RA        <NA>
##        ...      ...                  ...    ... .      ...      ...       ...       ...         ...         ...           ...               ...         ...
##   [160890]       2L     [759231, 759549]      + |  FlyBase      CDS      <NA>         0 FBgn0031275   GABA-B-R3   FBtr0343787      GABA-B-R3-RG        <NA>
##   [160891]       2L     [759679, 760910]      + |  FlyBase      CDS      <NA>         2 FBgn0031275   GABA-B-R3   FBtr0343787      GABA-B-R3-RG        <NA>
##   [160892]       2L     [773595, 773915]      + |  FlyBase      CDS      <NA>         0 FBgn0031276     CG12506   FBtr0078003        CG12506-RA        <NA>
##   [160893]       2L     [776542, 776838]      + |  FlyBase      CDS      <NA>         0 FBgn0040725     CG13946   FBtr0078004        CG13946-RA        <NA>
##   [160894]       2L     [779204, 779560]      + |  FlyBase      CDS      <NA>         0 FBgn0031277     CG13947   FBtr0078005        CG13947-RA        <NA>
##   -------
##   seqinfo: 25 sequences from an unspecified genome; no seqlengths
```



```r
# extract by gene symbol

fly_rRNA <- fly_gtf[grep("rRNA", fly_gtf$gene_symbol)]
```


```
## GRanges object with 416 ranges and 9 metadata columns:
##         seqnames               ranges strand |   source       type     score     phase     gene_id          gene_symbol transcript_id       transcript_symbol           #
##            <Rle>            <IRanges>  <Rle> | <factor>   <factor> <numeric> <integer> <character>          <character>   <character>             <character> <character>
##     [1]       2R [19766144, 19766278]      - |  FlyBase       gene      <NA>      <NA> FBgn0053353       5SrRNA:CR33353          <NA>                    <NA>        <NA>
##     [2]       2R [19766144, 19766278]      - |  FlyBase       rRNA      <NA>      <NA> FBgn0053353       5SrRNA:CR33353   FBtr0086345       5SrRNA:CR33353-RA        <NA>
##     [3]       2R [19766144, 19766278]      - |  FlyBase       exon      <NA>      <NA> FBgn0053353       5SrRNA:CR33353   FBtr0086345       5SrRNA:CR33353-RA        <NA>
##     [4]       2R [19765775, 19765909]      - |  FlyBase       gene      <NA>      <NA> FBgn0053354       5SrRNA:CR33354          <NA>                    <NA>        <NA>
##     [5]       2R [19765775, 19765909]      - |  FlyBase       rRNA      <NA>      <NA> FBgn0053354       5SrRNA:CR33354   FBtr0086346       5SrRNA:CR33354-RA        <NA>
##     ...      ...                  ...    ... .      ...        ...       ...       ...         ...                  ...           ...                     ...         ...
##   [412]        X [23291549, 23291671]      + |  FlyBase pseudogene      <NA>      <NA> FBgn0267523 5.8SrRNA-Psi:CR45863   FBtr0346897 5.8SrRNA-Psi:CR45863-RA        <NA>
##   [413]        X [23291549, 23291671]      + |  FlyBase       exon      <NA>      <NA> FBgn0267523 5.8SrRNA-Psi:CR45863   FBtr0346897 5.8SrRNA-Psi:CR45863-RA        <NA>
##   [414]        X [23291700, 23291729]      + |  FlyBase       gene      <NA>      <NA> FBgn0267524       2SrRNA:CR45864          <NA>                    <NA>        <NA>
##   [415]        X [23291700, 23291729]      + |  FlyBase       rRNA      <NA>      <NA> FBgn0267524       2SrRNA:CR45864   FBtr0346898       2SrRNA:CR45864-RA        <NA>
##   [416]        X [23291700, 23291729]      + |  FlyBase       exon      <NA>      <NA> FBgn0267524       2SrRNA:CR45864   FBtr0346898       2SrRNA:CR45864-RA        <NA>
##   -------
##   seqinfo: 25 sequences from an unspecified genome; no seqlengths
```
