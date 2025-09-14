/*
 * Trim adapters from reads
 */
process trim_using_cutadapt {

   tag "${id}:q > ${trim_qual}:l > ${min_length}"
   label "big_cpu" 
   
   // errorStrategy 'retry'
   // maxRetries 1

   publishDir( 
      "${params.outputs}/trimmed", 
      mode: 'copy',
   )

   input:
   tuple val( id ), path( reads, stageAs: "?/*" ), val( adapters5 ), val( adapters3 )
   tuple val( trim_qual ), val( min_length )

   output:
   tuple val( id ), path( "*.with-adapters.fastq.gz" ), emit: main
   tuple val( id ), path( "*.log" ), emit: logs

   script:
   def adapter_3prime_R = (adapters3[1].length() > 0 ? "-A '${adapters3[1]}'" : "")
   def adapter_5prime_R = (adapters5[1].length() > 0 ? "-G '${adapters5[1]}'" : "")
   def adapter_5prime_F = (adapters5[0].length() > 0 ? "-g '${adapters5[0]}'" : "")
   """
   for i in \$(seq 1 2)
   do
      cat */*_R"\$i"*.fastq.gz > ${id}_3p_R"\$i".fastq.gz
   done

   cutadapt \
		-a "${adapters3[0]}" \
      ${adapter_3prime_R} \
      --nextseq-trim ${trim_qual} \
      --minimum-length ${min_length} \
		--report full \
      --action trim \
      -j ${task.cpus} \
		-o ${id}_5p_R1.fastq.gz \
      -p ${id}_5p_R2.fastq.gz \
		${id}_3p_R?.fastq.gz \
   > ${id}.3p.cutadapt.log

   ADAPT5_ALL="${adapter_5prime_F}${adapter_5prime_R}"
   ADAPT5_LEN=\${#ADAPT5_ALL}
   if [ \$ADAPT5_LEN -gt 0 ]
   then
      cutadapt \
         ${adapter_5prime_F} ${adapter_5prime_R} \
         --report full \
         --action retain \
         --discard-untrimmed \
         --minimum-length ${min_length} \
         -j ${task.cpus} \
         -o ${id}_R1.with-adapters.fastq.gz \
         -p ${id}_R2.with-adapters.fastq.gz \
         ${id}_5p_R?.fastq.gz > ${id}.5p.cutadapt.log
   else
      for i in \$(seq 1 2)
      do
         mv ${id}_5p_R\$i.fastq.gz ${id}_R\$i.with-adapters.fastq.gz
      done
   fi

   rm ${id}_{3p,5p}_R?.fastq.gz
   """
}


process trim_nanopore_using_cutadapt {

   tag "${id}:q > ${trim_qual}:l > ${min_length}"
   label "big_mem" 

   publishDir( 
      "${params.outputs}/trimmed",
      mode: 'copy',
   )

   input:
   tuple val( id ), path( reads ), val( adapter5 ), val( adapter3 )
   tuple val( trim_qual ), val( min_length )

   output:
   tuple val( id ), path( "*.with-adapters.fastq.gz" ), emit: main
   tuple val( id ), path( "*.log" ), emit: logs

   script:
   """
   cat ${reads} > "${id}_5prime.fastq.gz"
   ERROR_RATE=0.35
   MIN_OVERLAP=15

   # `--revcomp` detects the strand
   cutadapt \
		-g '${adapter5}' \
      --error-rate \$ERROR_RATE \
      --overlap \$MIN_OVERLAP \
		-q ${trim_qual.toString().split(',')[0]} \
      --minimum-length ${min_length.toString().split(':')[0]} \
      --revcomp \
      -j ${task.cpus} \
		--report full \
      --action retain \
		-o "${id}_3prime.fastq.gz" \
		"${id}_5prime.fastq.gz" \
   > "${id}_5prime.cutadapt.log"

   cutadapt \
		-a '${adapter3}' \
      --error-rate \$ERROR_RATE \
      --overlap \$MIN_OVERLAP \
		--minimum-length ${min_length.toString().split(':')[0]} \
      -j ${task.cpus} \
		--report full \
      --action retain \
		-o "${id}.with-adapters.fastq.gz" \
		"${id}_3prime.fastq.gz" \
   > "${id}_3prime.cutadapt.log"

   """
}