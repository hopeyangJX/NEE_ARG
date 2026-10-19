# About
### This is the main script used for metagenomic data processing and statistical analysis in the following study:
Soil Acidification Mediates Nitrogen‑Driven Resistome Expansion and Offers a Validated Mitigation Target. *Unpublished*.

For simplicity, we take sample S1 as an example:
# 1 Metagenomic data preprocessing
## 1.1 quality control
We used fastp (https://github.com/opengene/fastp) for QC.
```bash
fastp \
    -i ${sample}_1.fastq.gz \
    -I ${sample}_2.fastq.gz \
    -o /path/to/clean/${sample}_1.clean.fastq.gz \
    -O /path/to/clean/${sample}_2.clean.fastq.gz \
    -q 30 \
    -t 20 \
    -h ${sample}_fastp.html \
    -j ${sample}_fastp.json
```

## 1.2 Assembly
We used MEGAHIT（https://github.com/voutcn/megahit) to assemble the contigs.

```bash
megahit \
    -t 40 \
    -m 0.3 \
    -1 ${sample}_1.clean.fastq.gz \
    -2 ${sample}_2.clean.fastq.gz \
    --k-list 21,29,39,59,79,99,119,141 \
    -o /path/to/megahit/${sample}.out \
    --out-prefix ${sample}
```

## 1.3 Prodigal
We used Prodigal （https://github.com/hyattpd/prodigal) for the prediction of Open reading frameworks (ORFs).
```bash
prodigal -p meta -q -f gff \
    -i ${sample}.contigs.fa \
    -o ${sample}.gff \
    -a ${sample}.faa \
    -d ${sample}.fna
```
## 1.4 geNomad
We use geNomad (https://github.com/apcamargo/genomad) to extract plasmid or viral contigs from assembled contigs for each treatment.
```bash
for fa in ${CONTIGS_DIR}/*.contig.fa; do
    sample=$(basename "${fa}" .contig.fa)
    nohup genomad end-to-end \
        --cleanup \
        --splits 8 \
        -t 20 \
        "${fa}" \
        ${GENOMAD_OUT_DIR}/${sample} \
        ${GENOMAD_DB} \
        > ${GENOMAD_OUT_DIR}/${sample}.log 2>&1 &
done
```

## 1.5 Binning
MetaWRAP（https://github.com/bxlab/metaWRAP) was employed for binning.
```bash
sample=<sample>
conda activate ${CONDA_ENV}
nohup metawrap binning \
    -o ${BIN_DIR}/${sample}.bin \
    -a ${CONTIGS_DIR}/${sample}.contigs.fa \
    --metabat2 \
    --maxbin2 \
    --concoct \
    ${READS_DIR}/${sample}_1.clean.fastq \
    ${READS_DIR}/${sample}_2.clean.fastq \
    -t ${THREADS} \
    > ${BIN_DIR}/${sample}.bin.log 2>&1 &
```

## 1.6 MAG dereplication
MAGs were dereplicated using dRep with 95% ANI and 60% completeness.
```bash
conda activate dRep
mkdir -p ${DREP_OUT_DIR}
nohup dRep dereplicate ${DREP_OUT_DIR} \
    -g ${REFINED_BIN_DIR}/*.fa \
    -p ${THREADS} \
    --P_ani 0.965 \
    --S_ani 0.965 \
    -strW 101 \
    --completeness 60 \
    > ${DREP_OUT_DIR}/drep.log 2>&1 &
```  

## 1.7 Genome annotation
We performed taxonomic classification of dereplicated MAGs using GTDB-Tk (https://github.com/Ecogenomics/GTDBTk).
```bash
nohup gtdbtk classify_wf \
    --genome_dir ${GENOME_DIR} \
    --out_dir ${GTDBTK_OUT_DIR} \
    --cpu 60 \
    -x fa \
    > ${GTDBTK_OUT_DIR}/gtdbtk.log 2>&1 &
```

# 2 ARG annotation
ARGs were annotated from predicted ORFs using CARD RGI and Resfams. The workflow includes CARD RGI annotation, extraction of RGI-annotated ORF IDs, Resfams annotation of unannotated ORFs, removal of efflux genes, best hit selection, merging CARD and Resfams annotations, mapping ARG ORFs to contig IDs, ORF coverage calculation, ARG abundance matrix construction, and normalization by sequencing depth.
## 2.1 CARD RGI annotation
We used RGI bwt (https://github.com/arpcard/rgi/blob/master/docs/rgi_bwt.rst) to align assembled contigs or MAGs against CARD’s canonical reference database.
```bash
conda activate rgi_6.0.8
mkdir -p /path/to/project/outfile
for faa in /path/to/project/faa/*.faa; do
    sample=$(basename "$faa" .faa)
    nohup rgi main --input_sequence "$faa" --output_file /path/to/project/outfile/${sample}.card_orf --local --clean --low_quality -t protein -a DIAMOND -n 8 &
done
```
## 2.2 Extract RGI-annotated ORF IDs
```bash
mkdir -p /path/to/project/outfile/orf_id
for f in /path/to/project/outfile/*.card_orf.txt; do
    awk 'BEGIN{FS="#"} NR>1{print $1}' "$f" > /path/to/project/outfile/orf_id/$(basename "$f")
done
```
## 2.3 Resfams annotation for unannotated ORFs
Filter unannotated ORFs
```bash
mkdir -p /path/to/project/unannotated_faa
for id_file in /path/to/project/outfile/orf_id/*.card_orf.txt; do
    sample=$(basename "$id_file" .card_orf.txt)
    awk 'NR==FNR{del[$0]=1;next} /^>/{split($0,p," #"); keep=!(substr(p[1],2) in del)} keep' "$id_file" /path/to/project/faa/${sample}.faa > /path/to/project/unannotated_faa/${sample}.filtered.faa
done
```
The remaining unannotated genes were filtered and subsequently annotated with Resfams core database (https://www.dantaslab.org/resfams)
```bash
conda activate HMMER
mkdir -p /path/to/project/resfams
for f in /path/to/project/unannotated_faa/*.filtered.faa; do
    s=$(basename "$f" .filtered.faa)
    hmmsearch --cut_ga --cpu 40 --tblout /path/to/project/resfams/$s.resfams /path/to/Resfams.hmm "$f"
done
```
## 2.4 Remove efflux genes and keep best hit
```bash
mkdir -p /path/to/project/de_resfams/best
for f in /path/to/project/resfams/*.resfams; do
    base=$(basename "$f" .resfams)
    awk '$3 !~ /efflux/' "$f" > /path/to/project/de_resfams/$base.resfams
    awk '/^#/{print > out} !/^#/{orf=$1; score=$6; if(!(orf in best)||score>best[orf]){best[orf]=score; line[orf]=$0}} END{for(orf in line) print line[orf] >> out}' out=/path/to/project/de_resfams/best/$base.resfams /path/to/project/de_resfams/$base.resfams
done
```

## 2.5 Merge CARD and Resfams annotations
 Extract Resfams and ARO IDs
```bash
mkdir -p /path/to/project/outfile/resfam_id /path/to/project/outfile/ARO_id
for f in /path/to/project/de_resfams/best/*.resfams; do
    base=$(basename "$f" .resfams)
    awk 'BEGIN{FS="[ ]+";OFS="\t"} /^#/{next} {print $1,$4}' "$f" > /path/to/project/outfile/resfam_id/$base.resfam_id.tsv
done
for f in /path/to/project/outfile/*.card_orf.txt; do
    base=$(basename "$f" .card_orf.txt)
    awk 'BEGIN{FS="\t";OFS="\t"} NR==1{print "ORF_ID","ARO";next} {split($1,a," # "); print a[1],$11}' "$f" > /path/to/project/outfile/ARO_id/$base.ARO.tsv
done
```

 Merge
```bash
mkdir -p /path/to/project/outfile/two_database
for res in /path/to/project/outfile/resfam_id/*.resfam_id.tsv; do
    base=$(basename "$res" .resfam_id.tsv)
    aro=/path/to/project/outfile/ARO_id/$base.ARO.tsv
    [ -f "$aro" ] || continue
    awk 'BEGIN{OFS="\t";print "ORF_ID","ID","database"} NR==FNR{print $1,$2,"resfams";next} FNR>1{print $1,$2,"ARO"}' "$res" "$aro" > /path/to/project/outfile/two_database/$base.tsv
done
```

## 2.6 Map ARG ORFs to contig IDs
```bash
mkdir -p /path/to/project/AMCC/arg
for f in /path/to/project/outfile/two_database/*.tsv; do
    b=$(basename "$f" .tsv)
    awk 'BEGIN{FS=OFS="\t"} NR==1{print $1,$2,$3,"contigID";next} {tmp=$1; sub(/_[0-9]+$/,"",tmp); print $1,$2,$3,tmp}' "$f" > /path/to/project/AMCC/arg/$b.tsv
done
```
## 2.7 Calculate ORF coverage with CoverM
```bash
conda activate coverM
mkdir -p /path/to/project/orf_bam /path/to/project/coverm
for fna in /path/to/project/fna/*.fna; do
    sample=$(basename "$fna" .fna)
    nohup coverm contig -1 /path/to/project/clean/${sample}_1.clean.fastq.gz -2 /path/to/project/clean/${sample}_2.clean.fastq.gz --reference "$fna" -t 40 --bam-file-cache-directory /path/to/project/orf_bam -o /path/to/project/coverm/$sample.fna.coverm.tsv &
done
```
## 2.8 Attach coverage to ARG ORFs and build matrix
```bash
mkdir -p /path/to/project/AMCC/arg/abundance /path/to/project/AMCC/arg/matrix
for f in /path/to/project/AMCC/arg/*.tsv; do
    fname=$(basename "$f" .tsv)
    cov=/path/to/project/coverm/$fname.fna.coverm.tsv
    [ -f "$cov" ] || continue
    awk -F'\t' 'BEGIN{OFS="\t"} NR==FNR{gsub(/\r/,"",$1);cov[$1]=$2;next} FNR==1{print $0,"Coverage_Mean";next} {gsub(/\r/,"",$1);print $0,($1 in cov?cov[$1]:0)}' "$cov" "$f" > /path/to/project/AMCC/arg/abundance/$fname.abundance.tsv
    awk -F'\t' 'BEGIN{OFS="\t"} FNR>1{sum[$2]+=$5} END{for(id in sum) print id,sum[id]}' /path/to/project/AMCC/arg/abundance/$fname.abundance.tsv > /path/to/project/AMCC/arg/matrix/$fname.matrix.tsv
done
```

## 2.9 Combine all samples and normalize
```bash
awk -F'\t' 'BEGIN{OFS=","} FNR==1{file_idx++;fname=FILENAME;sub(/.*\//,"",fname);sub(/\..*/,"",fname);headers[file_idx]=fname} {ids[$1]=1;data[$1,file_idx]=$2} END{printf "ID";for(i=1;i<=file_idx;i++)printf ",%s",headers[i];print "";for(id in ids){printf "%s",id;for(i=1;i<=file_idx;i++){val=data[id,i];printf ",%s",(val==""?0:val)}print ""}}' /path/to/project/AMCC/arg/matrix/*.tsv | sort -t',' -k1,1 > /path/to/project/AMCC/arg/matrix/final_merged_matrix.csv
```

# After transposing and adding sequencing depth in column 2:
```bash
awk -F',' 'BEGIN{OFS=","} NR==1{print;next} {for(i=3;i<=NF;i++) if($2!=0) $i=$i/$2; print}' /path/to/project/AMCC/arg/matrix/T_final_merged_matrix.csv > /path/to/project/AMCC/arg/matrix/final_merged_matrix.normalized.csv
```

# 3 MRG, ABG, and ROS annotation
MRGs, ABGs, and ROS genes were annotated from predicted ORFs using DIAMOND blastp against BacMet, MIBiG, and ROSBase, respectively. 
The common workflow is identical; only the reference database, filtering thresholds, and downstream summarization differ. CoverM (https://github.com/wwood/CoverM) was used to calculate ORF coverage.
## 3.1 DIAMOND blastp
```bash
P=/path/to/project; F=$P/faa; DB=/path/to/<tag>.dmnd; O=$P/<tag>_diamond; mkdir -p $O
for f in $F/*.faa; do s=$(basename $f .faa); nohup diamond blastp -d $DB -q $f -f 6 qseqid sseqid pident length qlen slen qcovhsp scovhsp evalue bitscore --very-sensitive --evalue 1e-5 --max-target-seqs 1 -o $O/${s}_raw.tsv > $O/${s}.log 2>&1 & done
```
## 3.2 Add header to DIAMOND output
```bash
P=/path/to/project; I=$P/<tag>_diamond; O=$I/<tag>; mkdir -p $O
for r in $I/*_raw.tsv; do s=$(basename $r _raw.tsv); { echo -e "qseqid\tsseqid\tpident\tlength\tqlen\tslen\tqcovhsp\tscovhsp\tevalue\tbitscore"; cat $r; } > $O/${s}_<tag>_diamond.tsv; done
```
## 3.3 Filter hits
Set PID, QCOV, and EV to gene-specific cutoffs for MRG (pident ≥ 60, qcovhsp ≥ 70, evalue ≤ 1e-5), ABG (pident ≥ 80, qcovhsp ≥ 70, evalue ≤ 1e-7), and ROS (pident ≥ 50, qcovhsp ≥ 70, evalue ≤ 1e-5), keeping one best hit per ORF.
```bash
P=/path/to/project; I=$P/<tag>_diamond/<tag>; O=$P/<tag>_diamond/list1; mkdir -p $O
for f in $I/*_<tag>_diamond.tsv; do s=$(basename $f _<tag>_diamond.tsv); awk -F'\t' -v pid=$PID -v qcov=$QCOV -v ev=$EV 'BEGIN{OFS="\t";print "Name","ResGeneID"} NR>1 && $3>=pid && $7>=qcov && $9<=ev && !seen[$1]{print $1,$2;seen[$1]=1}' $f > $O/${s}_<tag>_list1.tsv; done
```

## 3.4 Calculate ORF coverage with CoverM
```
coverm contig -1 ${sample}_1.clean.fastq.gz -2 ${sample}_2.clean.fastq.gz --reference ${sample}.fna -t 20 --bam-file-cache-directory orf_bam -o ${sample}.fna.coverm.tsv
```
## 3.5 Merge ORF coverage with annotation list
```bash
P=/path/to/project; L=$P/<tag>_diamond/list1; C=$P/coverm; O=$P/<tag>_diamond/coverm; mkdir -p $O
for f in $L/*_<tag>_list1.tsv; do s=$(basename $f _<tag>_list1.tsv); awk 'BEGIN{OFS="\t"} NR==FNR{c[$1]=$2;next}{print $0,(FNR==1?"Coverage_Mean":($1 in c?c[$1]:0))}' $C/${s}.fna.coverm.tsv $f > $O/${s}_<tag>_coverm.tsv; done
```

## 3.6 Summarize total coverage per sample
```bash
D=/path/to/project/<tag>_diamond/coverm; echo "Sample,Total_Coverage" > $D/<tag>_coverage_summary.csv
for f in $D/*_<tag>_coverm.tsv; do s=$(basename $f | cut -d_ -f1); awk -v s=$s 'NR>1{t+=$3}END{printf "%s,%.6f\n",s,t}' $f >> $D/<tag>_coverage_summary.csv; done
```


# 4 Co-localization of ARGs and MRGs
ARG–MRG co-localization was assessed on assembled contigs, including plasmid contigs identified by geNomad. 
The same workflow applies to both plasmid contigs and the full set of assembled contigs; only the input contig set differs. 
ARG and MRG ORFs were first mapped back to their contig IDs, and only contigs carrying both an ARG and an MRG were retained as ARG–MRG co-localized contigs (AMCC). 
Physical distances between the ARG and MRG on the same contig were then calculated, and co-localized pairs within 5 kb were used for downstream analysis.
## 4.1 Map ARG ORFs to contig IDs
```bash
P=/path/to/project; I=$P/genomad/plasmid/two_database; O=$P/genomad/plasmid/AMCC/arg; mkdir -p $O
for f in $I/*.tsv; do b=$(basename $f .tsv); awk 'BEGIN{FS=OFS="\t"} NR==1{print $1,$2,$3,"contigID";next}{tmp=$1; sub(/_[0-9]+$/,"",tmp); print $1,$2,$3,tmp}' $f > $O/$b.tsv; done
```
## 4.2 Map MRG ORFs to contig IDs
```bash
P=/path/to/project; I=$P/genomad/plasmid/bacmet_diamond/list1; O=$P/genomad/plasmid/AMCC/mrg; mkdir -p $O
for f in $I/*.tsv; do b=$(basename $f .tsv); awk 'BEGIN{FS=OFS="\t"} NR==1{print $1,$2,"contigID";next}{o=$1; sub(/_[^_]+$/,"",o); print $1,$2,o}' $f > $O/$b.tsv; done
```
## 4.3 Merge ARG and MRG by shared contig ID
```bash
P=/path/to/project; A=$P/genomad/plasmid/AMCC/arg; M=$P/genomad/plasmid/AMCC/mrg; O=$P/genomad/plasmid/AMCC/arg_mrg; mkdir -p $O
for a in $A/*.tsv; do s=$(basename $a .tsv); m=$M/${s}_bacmet_list1.tsv; [ -f $m ] || continue
  awk 'BEGIN{FS=OFS="\t"} NR==FNR{if(FNR>1){n=++c[$3]; N[$3","n]=$1; R[$3","n]=$2} next} FNR==1{print "ORF_ID","ID","database","contigID","Name","ResGeneID"; next} $4 in c{for(i=1;i<=c[$4];i++) print $1,$2,$3,$4,N[$4","i],R[$4","i]}' $m $a > $O/${s}_merged.tsv
done
```
## 4.4 Add contig length and filter by size
Contig length was extracted from the assembled contigs and used to retain AMCCs between 1 kb and 5 kb.
```bash
P=/path/to/project; M=$P/genomad/plasmid/AMCC/arg_mrg; F=$P/contigs; O=$P/genomad/plasmid/AMCC/sequence_depth; mkdir -p $O
for f in $M/*_merged.tsv; do s=$(basename $f _merged.tsv); fa=$F/${s}.contig.fna; [ -f $fa ] || continue
  awk -v fa=$fa 'BEGIN{FS=" ";OFS="\t"} FILENAME==fa{if(/^>/){if(n!="")L[n]=l; n=substr($1,2); l=0}else l+=length($0); next} ENDFILE&&FILENAME==fa{if(n!="")L[n]=l; next} FNR==1{print $0,"sequence_depth"; next} {print $0,($4 in L?L[$4]:"NA")}' $fa $f > $O/${s}_depth.tsv
done
P=/path/to/project; I=$P/genomad/plasmid/AMCC/sequence_depth; O=$P/genomad/plasmid/AMCC/5kb; mkdir -p $O
for f in $I/*_depth.tsv; do s=$(basename $f _depth.tsv); awk 'BEGIN{FS=OFS="\t"} NR==1{print;next} $7>=1000 && $7<=5000{print}' $f > $O/${s}_5kb.tsv; done
```
## 4.5 Extract ARG and MRG physical coordinates
Physical coordinates of ARG and MRG ORFs were extracted from CARD RGI, Resfams, and BacMet outputs, then merged per ORF.
# ARG coordinates from CARD RGI
```bash
P=/path/to/project; I=$P/genomad/plasmid/outfile; O=$P/genomad/plasmid/outfile/ARO_ID_allname; mkdir -p $O
for f in $I/*.card_orf.txt; do s=$(basename $f .card_orf.txt); awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "ORF_ID","ARO","Start","End","Strand"; next} {split($1,p," # "); split($1,i," "); print i[1],$11,p[2],p[3],p[4]}' $f > $O/${s}_ARO_with_coords.txt; done
```
# ARG coordinates from Resfams
```bash
P=/path/to/project; I=$P/genomad/plasmid/de_resfams/best; O=$P/genomad/plasmid/outfile/resfam_id_allname; mkdir -p $O
for f in $I/*.resfams; do s=$(basename $f .resfams); awk 'BEGIN{OFS="\t"; print "ORF_ID","Resfam","Start","End","Strand"} /^#/{next} {split($0,h," # "); print $1,$4,h[2],h[3],h[4]}' $f > $O/${s}_resfams_with_coords.txt; done
```
# MRG coordinates from BacMet
```bash
P=/path/to/project; I=$P/genomad/plasmid/bacmet_diamond/list1_allname; O=$P/genomad/plasmid/AMCC/mrg_PID; mkdir -p $O
for f in $I/*_bacmet_list1_with_coords.tsv; do s=$(basename $f _bacmet_list1_with_coords.tsv); awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print "ORF_ID","MRG_Hit","Contig","Start","End","Strand"; next} {print $1,$2,$(NF-3),$(NF-2),$(NF-1),$NF}' $f > $O/${s}_mrg_PID.tsv; done
```
## 4.6 Merge ARG and MRG coordinates with AMCC
# Merge ARG coordinates (CARD + Resfams)
```bash
P=/path/to/project; R=$P/genomad/plasmid/outfile/resfam_id_allname; A=$P/genomad/plasmid/outfile/ARO_ID_allname; O=$P/genomad/plasmid/AMCC/arg/two_pid; mkdir -p $O
for r in $R/*_resfams_with_coords.txt; do s=$(basename $r _resfams_with_coords.txt); a=$A/${s}_ARO_with_coords.txt; [ -f $a ] || continue
  { echo -e "ORF_ID\tARO_Resfam\tStart\tEnd\tStrand"; awk -F'\t' 'BEGIN{OFS="\t"} NR>1{print $1,$2,$3,$4,$5}' $r; awk -F'\t' 'BEGIN{OFS="\t"} NR>1{print $1,$2,$3,$4,$5}' $a; } > $O/${s}.tsv
done
```
# Attach ARG coordinates to ARG contig table
```bash
P=/path/to/project; B=$P/genomad/plasmid/AMCC/arg; PID=$P/genomad/plasmid/AMCC/arg/two_pid; O=$P/genomad/plasmid/AMCC/arg/combine_pid; mkdir -p $O
for f in $B/*.tsv; do s=$(basename $f .tsv); p=$PID/${s}.tsv; [ -f $p ] || continue
  awk -F'\t' 'BEGIN{OFS="\t"} FILENAME==ARGV[1]{if(FNR>1) d[$1]=$2"\t"$3"\t"$4"\t"$5; next} FNR==1{print "ARG_ID","ID","database","contigID","ARO_Resfam","ARG_Start","ARG_End","ARG_Strand"; next} $1 in d{print $1,$2,$3,$4,d[$1]}' $p $f > $O/${s}_combined_ARG.tsv
done
```
# Attach MRG coordinates to MRG contig table
```bash
P=/path/to/project; B=$P/genomad/plasmid/AMCC/mrg; PID=$P/genomad/plasmid/AMCC/mrg_PID; O=$P/genomad/plasmid/AMCC/mrg/MRG_ID; mkdir -p $O
for f in $B/*_bacmet_list1.tsv; do s=$(basename $f _bacmet_list1.tsv); p=$PID/${s}_mrg_PID.tsv; [ -f $p ] || continue
  awk -F'\t' 'BEGIN{OFS="\t"} FILENAME==ARGV[1]{if(FNR>1) d[$1]=$3"\t"$4"\t"$5"\t"$6; next} FNR==1{print "MRG","ResGeneID","contigID","MRG_Start","MRG_End","MRG_Strand"; next} $1 in d{print $1,$2,d[$1]}' $p $f > $O/${s}.tsv
done
```
## 4.7 Merge ARG and MRG on the same contig with coordinates
```bash
P=/path/to/project; M=$P/genomad/plasmid/AMCC/mrg/MRG_ID; A=$P/genomad/plasmid/AMCC/arg/combine_pid; O=$P/genomad/plasmid/AMCC/arg_mrg_pid; mkdir -p $O
for m in $M/*.tsv; do s=$(basename $m .tsv); a=$A/${s}_combined_ARG.tsv; [ -f $a ] || continue
  awk -F'\t' 'BEGIN{OFS="\t"} FILENAME==ARGV[1]{if(FNR>1) d[$4]=$0; next} FNR==1{print $0,"ARG_ID","ID","database","contigID_ARG","ARO_Resfam","ARG_Start","ARG_End","ARG_Strand"; next} $3 in d{print $0,d[$3]}' $a $m > $O/${s}_ARG_MRG_merged.tsv
done
```
## 4.8 Calculate physical distance between ARG and MRG
Physical distance was defined as the gap between the two ORFs: max(start) − min(end) if they do not overlap, and 0 if they overlap.
```bash
P=/path/to/project; I=$P/genomad/plasmid/AMCC/arg_mrg_pid; O=$P/genomad/plasmid/AMCC/arg_mrg_distance; mkdir -p $O
for f in $I/*_ARG_MRG_merged.tsv; do s=$(basename $f _ARG_MRG_merged.tsv)
  awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print $0,"Distance_bp"; next} {ms=$4; me=$5; as=$12; ae=$13; mx=(ms>as?ms:as); mn=(me<ae?me:ae); print $0,(mx>mn?mx-mn:0)}' $f > $O/${s}_ARG_MRG_with_dist.tsv
done
```
## 4.9 Keep co-localized pairs within 5 kb
```bash
P=/path/to/project; I=$P/genomad/plasmid/AMCC/arg_mrg_distance; O=$P/genomad/plasmid/AMCC/arg_mrg_distance/5kb; mkdir -p $O
for f in $I/*_ARG_MRG_with_dist.tsv; do awk -F'\t' 'BEGIN{OFS="\t"} NR==1{print;next} $NF<=5000{print}' $f > $O/$(basename $f); done
```
## 4.10 Calculate AMCC coverage
Contig coverage was attached to each co-localized ARG–MRG pair, and unique contigs were summed to obtain the total AMCC coverage per sample.
```
P=/path/to/project; I=$P/genomad/plasmid/AMCC/arg_mrg_distance/5kb; C=$P/coverm_contigs; O=$I/cv; mkdir -p $O
for f in $I/*_ARG_MRG_with_dist.tsv; do s=$(basename $f _ARG_MRG_with_dist.tsv); c=$C/${s}.contig.coverm.tsv; [ -f $c ] || continue
  awk -F'\t' 'BEGIN{OFS="\t"} FILENAME==ARGV[1]{if(FNR>1) v[$1]=$2; next} FNR==1{print $0,"Contig_Mean_Cov"; next} {print $0,($3 in v?v[$3]:"NA")}' $c $f > $O/${s}_ARG_MRG_with_cv.tsv
done
echo "Sample,Total_Contig_Mean_Cov" > $O/Contig_Mean_Cov_5KB.csv
for f in $O/*_ARG_MRG_with_cv.tsv; do s=$(basename $f | awk -F'_' '{print $1"_"$2}'); awk -F'\t' -v s=$s 'NR>1{if(!seen[$3]++){sum+=$NF}} END{printf "%s,%.6f\n",s,sum}' $f >> $O/Contig_Mean_Cov_5KB.csv; done
```
## 4.11 Build AMCC heatmap matrix
Each unique ARG–MRG pair was used as a row, and each sample as a column. Coverage values were summed when duplicated, and zeros were filled for missing pairs.
```bash
python3 -c "
import glob, os
D='/path/to/project/genomad/plasmid/AMCC/arg_mrg_distance/5kb/heatmap/coverm'
O=D+'/ARG_MRG_CV_matrix.csv'
fs=sorted(glob.glob(os.path.join(D,'*.tsv')))
hs=[os.path.basename(f).split('.')[0] for f in fs]
S=set(); d={}
for i,f in enumerate(fs):
    for j,l in enumerate(open(f)):
        if j==0: continue
        p=l.strip().split('\t')
        if len(p)>=8:
            S.add(p[6]); d[(p[6],i)]=d.get((p[6],i),0)+float(p[7])
open(O,'w').write('ARG_MRG_ID,'+','.join(hs)+'\n'+'\n'.join(a+','+','.join(str(d.get((a,i),0)) for i in range(len(fs))) for a in sorted(S)))
"
```
## 4.12 Normalize AMCC matrix by sequencing depth
After transposing the matrix and adding sample sequencing depth in column 2, each sample column was divided by its sequencing depth.
```bash
python3 -c "
import csv
F='/path/to/project/genomad/plasmid/AMCC/arg_mrg_distance/5kb/heatmap/coverm/T_AMCC_all_merged.csv'
O='/path/to/project/genomad/plasmid/AMCC/arg_mrg_distance/5kb/heatmap/coverm/T_P_5KBAMCC_normalized.csv'
with open(F) as fi, open(O,'w',newline='') as fo:
    w=csv.writer(fo)
    for i,r in enumerate(csv.reader(fi)):
        if i==0: w.writerow(r); continue
        try: g=float(r[1])
        except: g=0.0
        w.writerow([r[0],r[1]]+['0']*(len(r)-2) if g==0 else [r[0],r[1]]+[str(float(v)/g) if v else '0' for v in r[2:]])
"
```
Note: This workflow applies to both plasmid contigs and the full set of assembled contigs. To run it on assembled contigs, replace genomad/plasmid with the corresponding directory containing the assembled contigs and their annotation outputs.
# 5 Human pathogen annotation
Human pathogen–associated contigs were identified by taxonomic classification with Kraken2 (https://github.com/DerrickWood/kraken2), filtered against a human pathogen list, and quantified using CoverM.
The workflow includes contig taxonomy, coverage calculation, pathogen filtering, species-level coverage summation, and matrix construction.
## 5.1 Taxonomic classification of contigs with Kraken2
```bash
conda activate kraken2_2.1.6
mkdir -p /path/to/project/all_kraken2
cd /path/to/project/contigs
nohup bash -c 'for fa in *.fa; do sample="${fa%.fa}"; kraken2 --db /path/to/kraken2_DB_core_NT "$fa" --threads 60 --use-names --output "/path/to/project/all_kraken2/${sample}.kraken2"; done' 2>/dev/null &
```
## 5.2 Calculate contig coverage with CoverM
```bash
conda activate coverM
mkdir -p /path/to/project/coverm_contigs
for fa in /path/to/project/contigs/*.contigs.fa; do
    sample=$(basename "$fa" .contigs.fa)
    r1="/path/to/project/clean/${sample}.clean.fq1.gz"
    r2="/path/to/project/clean/${sample}.clean.fq2.gz"
    out="/path/to/project/coverm_contigs/${sample}.contig.coverm.tsv"
    coverm contig -1 "$r1" -2 "$r2" --reference "$fa" -t 50 --bam-file-cache-directory /path/to/project/orf_bam -o "$out"
done
```
## 5.3 Keep classified contigs
```bash
mkdir -p /path/to/project/all_kraken2/classified
cd /path/to/project/all_kraken2
nohup bash -c 'for f in *.kraken2; do grep "^C" "$f" > "classified/${f%.kraken2}.classified.txt"; done' 2>/dev/null &
```
## 5.4 Extract species names
```bash
mkdir -p /path/to/project/all_kraken2/species
cd /path/to/project/all_kraken2/classified
for f in *.classified.txt; do
  awk -F'\t' -v OFS='\t' '{split($3,a," \\("); $(NF+1)=a[1]; print}' "$f" > ../species/"${f%.classified.txt}.species.txt"
done
```
## 5.5 Filter human pathogen–associated contigs
A human pathogen list (human_pathogen.txt) is required, with species name in column 1 and pathogen type in column 2. 
```bash
mkdir -p /path/to/project/all_kraken2/human_pathogens
IN_DIR="/path/to/project/all_kraken2/species"
PAT_FILE="/path/to/project/human_pathogen.txt"
OUT_DIR="/path/to/project/all_kraken2/human_pathogens"
cd "$IN_DIR"
for f in *.species.txt; do
  [ "$f" = "human_pathogen.txt" ] && continue
  awk -F'\t' -v OFS='\t' 'FNR==NR{pathogen[$1]=$2; next} $6 in pathogen{print $0, pathogen[$6]}' "$PAT_FILE" "$f" > "$OUT_DIR/$f"
done
```
## 5.6 Extract pathogen contig IDs
```bash
mkdir -p /path/to/project/all_kraken2/hp_pathogen_contig_id
nohup find /path/to/project/all_kraken2/human_pathogens -name "*.species.txt" -exec sh -c '
out="/path/to/project/all_kraken2/hp_pathogen_contig_id/$(basename "$1" .txt).id.txt"
awk "BEGIN{OFS=\"\t\"; print \"pathogens_contig_id\", \"pathogens\"} {print \$2, \$(NF-3) \"_\" \$(NF-2)}" "$1" > "$out"
' sh {} \; &
```
## 5.7 Calculate pathogen coverage
```
mkdir -p /path/to/project/all_kraken2/hp_coverage
for id_file in /path/to/project/all_kraken2/hp_pathogen_contig_id/*.id.txt; do
    sample=$(basename "$id_file" .species.id.txt)
    coverm_file="/path/to/project/coverm_contigs/${sample}.contig.coverm.tsv"
    out_file="/path/to/project/all_kraken2/hp_coverage/${sample}.pathogen.coverage.tsv"
    awk 'NR==FNR{ids[$1]=$2; next} $1 in ids {print $0 "\t" ids[$1]}' "$id_file" "$coverm_file" > "$out_file"
done
```

## 5.8 Sum coverage by species
```bash
mkdir -p /path/to/project/all_kraken2/hp_species_coverage
for f in /path/to/project/all_kraken2/hp_coverage/*.tsv; do
    fname=$(basename "$f")
    awk -F'\t' 'BEGIN{OFS="\t"} {sum[$3]+=$2} END{for(i in sum) print i, sum[i]}' "$f" > "/path/to/project/all_kraken2/hp_species_coverage/${fname%.tsv}.merged.tsv"
done
```
## 5.9 Build species abundance matrix
```bash
python3 -c "
import glob, os
D='/path/to/project/all_kraken2/hp_species_coverage'
os.chdir(D)
fs=sorted(glob.glob('*.merged.tsv'))
hs=[f[:-11] for f in fs]
S=set(); d={}
for i,f in enumerate(fs):
    for l in open(f):
        p=l.strip().split('\t')
        if len(p)>=2 and p[0].lower() not in ('id','species',''):
            S.add(p[0]); d[(p[0],i)]=p[1]
open('final_merged_matrix.csv','w').write('ID,'+','.join(hs)+'\n'+'\n'.join(s+','+','.join(d.get((s,i),'0') for i in range(len(fs))) for s in sorted(S)))
"
```
## 5.10 Normalize by sequencing depth
After transposing the matrix and adding sample sequencing depth in column 2, the resulting matrix was further normalized by dividing each sample column by its corresponding sequencing depth.
```bash
awk -F',' 'BEGIN{OFS=","} NR==1{print;next} {for(i=3;i<=NF;i++) if($2!=0) $i=$i/$2; print}' /path/to/project/all_kraken2/hp_species_coverage/T_final_merged_matrix.csv > /path/to/project/all_kraken2/hp_species_coverage/hp_normalized.csv
```
