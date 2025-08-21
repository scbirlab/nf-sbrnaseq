process fetch_genome_from_NCBI {

   tag "${accession}"
   label 'some_mem'

   input:
   val accession

   output:
   tuple val( accession ), path( "all-nucleotides.fna" ), path( "all-annotations.gff" )

   script:
   """
   set -euox pipefail
   ACCESSIONS=\$(echo "${accession}" | tr '+' ' ')
   echo "\$ACCESSIONS"
   WEB_ROOT="https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession"
   WEB_TAIL="download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED"
   for acc in \$ACCESSIONS
   do
      curl -L "\$WEB_ROOT/\$acc/\$WEB_TAIL" \
      > \$acc"_genome-out"
      unzip -o \$acc"_genome-out" "ncbi_dataset/data/\$acc/"{"\$acc"_*_genomic.fna,*.gff}
      mv ncbi_dataset/data/*/"\$acc"_*_genomic.fna \$acc.fna
      mv ncbi_dataset/data/*/*.gff \$acc.gff
   done

   cat *.fna > all-nucleotides.fna
   cat *.gff > all-annotations.gff
   """
}


// Get FASTQ
process fetch_FASTQ_from_SRA {

   tag "${sample_id}-${sra_run_id}" 

   label 'big_mem'
   time '24 h'

   input:
   tuple val( sample_id ), val( sra_run_id )
   secret 'NCBI_API_KEY'

   output:
   tuple val( sample_id ), path( "*.with-idx_R?.fastq.gz" )

   script:
   """
   NCBI_API_KEY=\$NCBI_API_KEY \
   fastq-dump \
      --read-filter pass \
      --origfmt --defline-seq '@rd.\$si:\$sg:\$sn' \
      --split-3 ${sra_run_id}

   for i in \$(seq 1 2)
   do
      if [ \$i -eq 1 ]
      then
         for f in *_\$i.fastq
         do
            awk -F: 'NR%4==1 { a=\$2; alen=length(a); print \$0 } NR%4==2 { print a \$0 } NR%4==3 { print "+" } NR%4==0 { s = sprintf("%*s", alen, ""); print gensub(".", "F", "g", s) \$0 }' \
               \$f \
            | gzip -v --best \
            > \$(basename \$f .fastq).with-idx_R\$i.fastq.gz
         done
      else
         for f in *_\$i.fastq
         do
            awk -F: 'NR%4<3 { print \$0 } NR%4==3 { print "+" }' \
               \$f \
            | gzip -v --best \
            > \$(basename \$f .fastq).with-idx_R\$i.fastq.gz
         done
      fi
   done
   rm *.fastq
   """

   stub:
   """
   NCBI_API_KEY=\$NCBI_API_KEY \
   fastq-dump \
      -X 1000000 \
      --read-filter pass \
      --origfmt --defline-seq '@rd.\$si:\$sg:\$sn' \
      --split-3 ${sra_run_id}

   for i in \$(seq 1 2)
   do
      if [ \$i -eq 1 ]
      then
         for f in *_\$i.fastq
         do
            awk -F: 'NR%4==1 { a=\$2; alen=length(a); print \$0 } NR%4==2 { print a \$0 } NR%4==3 { print "+" } NR%4==0 { s = sprintf("%*s", alen, ""); print gensub(".", "F", "g", s) \$0 }' \
               \$f \
            | gzip -v --best \
            > \$(basename \$f .fastq).with-idx_R\$i.fastq.gz
         done
      else
         for f in *_\$i.fastq
         do
            awk -F: 'NR%4<3 { print \$0 } NR%4==3 { print "+" }' \
               \$f \
            | gzip -v --best \
            > \$(basename \$f .fastq).with-idx_R\$i.fastq.gz
         done
      fi
   done
   """
}