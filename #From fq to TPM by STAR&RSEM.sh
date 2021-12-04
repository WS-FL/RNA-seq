#Index build
STAR \
--runMode genomeGenerate \ 
--genomeDir ./ \ #Outdir
--runThreadN 10 \ #Threads
--genomeFastaFiles dir/hg38.fa \
--sjdbGTFfile dir/gencode.v38.annotation.gtf \
1>log.txt

#Use STAR to get bam
STAR \
--runThreadN 10 \
--runMode alignReads \
--quantMode TranscriptomeSAM GeneCounts \
--outSAMtype BAM SortedByCoordinate \ #It will be smaller than None
--outSAMunmapped None \
--genomeDir index/Prefix \
--readFilesIn ${dir}/${i}_1_val_1.fq \
              ${dir}/2_trim/${i}_2_val_2.fq \
--outFileNamePrefix $i 1>log.txt
#more choices
#--twopassMode Basic \ #先按索引进行第一次比对，而后把第一次比对发现的新剪切位点信息加入到索引中进行第二次比对。这个参数可以保证更精准的比对情况，但是费时也费内存。
#--readFilesCommand zcat \ #说明你的fastq文件是压缩形式的，就是.gz结尾的，不加的话会报错

#to get counts table
rsem-calculate-expression  --paired-end \
--no-bam-output --alignments -p 15 \ #Threads INT
-q SRRXXXXAligned.toTranscriptome.out.bam \ #from STAR toTranscriptome
../../rsem_test/Name_Prefix \ #index directory
rsemtest/Name_Prefix #Output Dir & Name Prefix '.genes.result' is Default

#merge the counts table 
rsem-generate-data-matrix *.result >output.matrix #only includes count matrix