process fetch_genome_from_NCBI {

   tag "${accession}"
   label 'some_mem'

   publishDir( 
      "${params.outputs}/genome", 
      mode: 'copy',
      saveAs: { "${accession}.${it}"},
   )

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
   WEB_TAIL="download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&hydrated=FULLY_HYDRATED"
   for acc in \$ACCESSIONS
   do
      curl -v "\$WEB_ROOT/\$acc/\$WEB_TAIL" \
         -o \$acc"_genome-out"
      unzip -o \$acc"_genome-out" "ncbi_dataset/data/\$acc/"{"\$acc"_*_genomic.fna,*.gff}
      mv ncbi_dataset/data/*/"\$acc"_*_genomic.fna \$acc.fna
      mv ncbi_dataset/data/*/*.gff \$acc.gff
   done

   cat *.fna > all-nucleotides.fna
   cat *.gff > all-annotations.gff
   """
}


// Get FASTQ
process prefetch_from_SRA {

   tag "${id}:${sra_run_id}" 

   input:
   tuple val( id ), val( sra_run_id )
   secret 'NCBI_API_KEY'

   output:
   tuple val( id ), path( "${sra_run_id}/" )

   script:
   """
   set -euxo pipefail 

   prefetch ${sra_run_id} --max-size u

   """

   stub:
   """
   mkdir "${sra_run_id}"

   """
}


// Get FASTQ
process download_FASTQ_from_SRA {

   tag "${id}:${sra_run_id}" 

   label 'big_cpu'

   input:
   tuple val( id ), path( sra_run_id )
   secret 'NCBI_API_KEY'

   output:
   tuple val( id ), path( "fastq/${sra_run_id}_{1,2}.fastq.gz", arity: 2 )

   script:
   """
   set -euxo pipefail
   fasterq-dump \
      --progress \
      --seq-defline '@\$ac:rd.\$si:\$sg:\$sn' \
      --qual-defline '+' \
      --threads ${task.cpus} \
      --outdir fastq \
      "${sra_run_id}"

   fastq_files=(fastq/*_{1,2}.fastq)
   for f in "\${fastq_files[@]}"
   do
      pigz -v -p ${task.cpus} "\$f"
   done

   """

   stub:
   """
   set -euxo pipefail 
   
   fastq-dump \
      -X 1000000 \
      --read-filter pass \
      --origfmt --defline-seq '@${sra_run_id}:rd.\$si:\$sg:\$sn' \
      --split-3 ${sra_run_id}

   fastq_files=(fastq/*_{1,2}.fastq)
   for f in "\${fastq_files[@]}"
   do
      pigz -v -p ${task.cpus} "\$f"
   done

   """
}

process prepend_reads_with_barcodes {

   tag "${id}" 

   publishDir( 
      "${params.outputs}/sra", 
      mode: 'copy',
      // saveAs: { "${id}.${it}" },
   )

   input:
   tuple val( id ), path( fastq )

   output:
   tuple val( id ), path( "*.with-idx_R{1,2}.fastq.gz", arity: 2 )

   script:
   """
   set -euxo pipefail

   for i in \$(seq 1 2)
   do
      if [ \$i -eq 1 ]
      then
         for f in *_\$i.fastq.gz
         do
            zcat "\$f" \
            | awk -F: '
               NR % 4 == 1 { 
                  a=\$3
                  alen=length(a)
                  print \$0
                  next
               } 
               NR % 4 == 2 { print a \$0; next }
               NR % 4 == 3 { print \$0; next } 
               NR % 4 == 0 { 
                  s = sprintf("%*s", alen, "")
                  gsub(/./, "F", s)
                  print s \$0 
               }
               ' \
            | pigz -v -p ${task.cpus} --stdout \
            > \$(basename "\$f" _\$i.fastq.gz).with-idx_R\$i.fastq.gz
         done
      else
         for f in *_\$i.fastq.gz
         do
            cp "\$f" \$(basename "\$f" _\$i.fastq.gz).with-idx_R\$i.fastq.gz
         done
      fi
   done

   """
}
