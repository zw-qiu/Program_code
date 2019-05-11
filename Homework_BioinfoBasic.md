Homework

1.2 Practice Guide

（1）第4列（$4）表示基因或转录本在参考序列上的起始位置，第5列（$5）表示基因或转录本在参考序列上的终止位置。Exon的长度应该是$5-$4+1。

（2）cat 1.gtf | awk '$1=="XI" && $3=="CDS"  {split ($10,x,";");name=x[1];gsub("\\"","",name);print $5,name}' | sort -n | tail -10 

| End_CDS of chr XI | gene_id |
| :---------------: | :-----: |
|      632798       | YKR097W |
|      635179       | YKR098C |
|      638283       | YKR099W |
|      639968       | YKR100C |
|      642501       | YKR101W |
|      649862       | YKR102W |
|      656733       | YKR103W |
|      657753       | YKR104W |
|      660464       | YKR105C |
|      663286       | YKR106W |

（3）grep -v '^#' 1.gtf | awk '$1="IV" {print $3,$2}' | sort |uniq -c 

|   .gtf_$3   | .gtf_$2 |
| :---------: | :-----: |
|     CDS     | ensembl |
|    exon     | ensembl |
|    gene     | ensembl |
| start_codon | ensembl |
| stop_codon  | ensembl |
| transcript  | ensembl |

