process fetch_UMIcollapse {

   input:
   val umicollapse_repo

   output:
   path "UMICollapse"

   script:
   """
   git clone ${umicollapse_repo} #https://github.com/siddharthab/UMICollapse.git
   UMI_LIB_PATH=UMICollapse/lib
   mkdir -p \$UMI_LIB_PATH
   curl -L https://repo1.maven.org/maven2/com/github/samtools/htsjdk/2.19.0/htsjdk-2.19.0.jar > \$UMI_LIB_PATH/htsjdk-2.19.0.jar
   curl -L https://repo1.maven.org/maven2/org/xerial/snappy/snappy-java/1.1.7.3/snappy-java-1.1.7.3.jar > \$UMI_LIB_PATH/snappy-java-1.1.7.3.jar
   """
}


process Prepend_Barcode {

   tag "${id}" 

   // label 'big_mem'
   // errorStrategy 'retry'
   // maxRetries 2

   // publishDir( 
   //    "${params.outputs}/umicollapse", 
   //    mode: 'copy',
   // )
   
   input:
   tuple val( id ), path( bamfile )
   val paired

   output:
   tuple val( id ), path( "prepended.bam" )

   script:
   """
   set -euox pipefail

   samtools view -h ${paired ? "" : "-f128"} "${bamfile[1]}" \
   | awk -F'\\t' -v OFS='\\t' '
      /^@/ { print; next }
      !/^@/ {
         delete a
         split(\$1, a, "_"); 
         \$1 = a[1] "_" a[2] "_" a[2] a[3]; 
         print \$0, "CB:Z:" a[2];
         next
      }
   ' | samtools view -h -bS - -o "prepended.bam"

   """
}


process Shard_for_UMIcollapse {

   tag "${id}" 

   label 'big_cpu'
   // errorStrategy 'retry'
   // maxRetries 2

   // publishDir( 
   //    "${params.outputs}/umicollapse", 
   //    mode: 'copy',
   // )
   
   input:
   tuple val( id ), path( bamfile )

   output:
   tuple val( id ), path( "shards/_header.sam" ), path( "shards/shard-*.bam" )

   script:
   """
   set -euox pipefail

   mkdir -p shards
   samtools view -H "${bamfile}" > shards/_header.sam
   #samtools sort "${bamfile}" \
   #   -@${task.cpus} -m ${Math.round(Math.floor(task.memory.getMega() * 0.8 / task.cpus))}M \
   #   -t CB \
   #   -o sorted.bam

   samtools view -@${task.cpus} "${bamfile}" \
   | awk -F'\\t' -v OFS='\\t' -v dir="shards" '
      {
         cb="";
         for(i=12;i<=NF;i++) if(\$i ~ /^CB:Z:/) { 
            split(\$i,a,":"); cb=a[3]; 
            break 
         }
         if(cb == "") next
         prefix = substr(cb, 1, 12); 
         if(prefix == "") prefix = "_no-prefix"
         print >> (dir "/shard-" prefix ".sam")
      }
   '

   for f in shards/shard-*.sam
   do
      name=\$(basename "\$f" .sam)
      cat shards/_header.sam "\$f" \
      | samtools view \
         -@${task.cpus} -bS \
         -o "shards/\${name}.bam" -
      rm "\$f"
   done

   """
}


process UMIcollapse {

   tag "${id}:${bamfile}" 

   label 'big_mem'
   // errorStrategy 'retry'
   // maxRetries 2

   publishDir( 
      "${params.outputs}/umicollapse/logs", 
      mode: 'copy',
      pattern: "*.log",
      saveAs: { "${id}.${bamfile.simpleName}.${it}" }
   )
   
   input:
   tuple val( id ), path( sam_header ), path( bamfile ), path( umicollapse_repo )
   val paired

   output:
   tuple val( id ), path( "umicollapse.bam" ), emit: main
   path "*.log", emit: logs

   script:
   """
   set -euox pipefail

   #samtools reheader "${sam_header}" "${bamfile}" > with-header.bam

   # sort by coordinate and stamp header as coordinate
   samtools sort \
      -@${task.cpus} \
      -m ${Math.round(Math.floor(task.memory.getMega() * 0.8 / task.cpus))}M \
      -o sorted.bam \
      "${bamfile}"
   # ensure @HD SO:coordinate (reheader if needed)
   samtools view -H sorted.bam \
   | awk '
      BEGIN { done = 0 } 
      /^@HD/ && \$0 ~ /SO:/ { sub(/SO:[^ \\t]+/,"SO:coordinate"); done=1 } 
      { print } 
      END { if(!done) print "@HD\\tVN:1.6\\tSO:coordinate" }
   ' \
   > header.sam
   samtools reheader header.sam sorted.bam > SOcoord.bam
   samtools index SOcoord.bam

   java -jar \
      -Xmx${Math.round(task.memory.getGiga() * 0.8)}G \
      -Xss1024m \
         "${umicollapse_repo}/umicollapse.jar" \
         bam ${paired ? "--paired" : ""} \
         --tag \
         -i SOcoord.bam \
         -o umicollapse.bam \
   >> "umicollapse.log" 2>&1
   rm SOcoord.bam header.sam sorted.bam

   if [ ! -e "umicollapse.log" ]
   then
      echo "There were no reads in ${bamfile}" > "umicollapse.log"
   fi

   """
}


process Concat_UMIcollapse {

   tag "${id}" 
   label 'big_cpu'
   // errorStrategy 'retry'
   // maxRetries 2

   publishDir( 
      "${params.outputs}/umicollapse/bam", 
      mode: 'copy',
      saveAs: { "${id}.${it}" }
   )
   
   input:
   tuple val( id ), path( original_bam ), path( bamfiles, stageAs: 'bamfiles-??????/dedup.bam' )

   output:
   tuple val( id ), path( "umicollapse.bam" ), emit: main

   script:
   """
   set -euox pipefail

   if ls *.bam > /dev/null 2>&1
   then
      samtools merge -@${task.cpus} \
         -h "${original_bam}" \
         -o merged.dedup.bam \
         bamfiles-*/dedup.bam
   else
      samtools view -H "${original_bam}" -o merged.dedup.bam
   fi

   samtools sort \
      -@${task.cpus} \
      -m ${Math.round(Math.floor(task.memory.getMega() * 0.8 / task.cpus))}M \
      -o "umicollapse.bam" \
      merged.dedup.bam

   samtools index "umicollapse.bam"
   rm merged.dedup.bam
   
   #start_count=\$(samtools view -c "${original_bam}")
   #end_count=\$(samtools view -c "umicollapse.bam")
   #if [ "\$start_count" -ne "\$end_count" ]
   #then
   #   >&2 echo "Lost some reads during duplicate tagging!"
   #   >&2 echo "- Initial count: \$start_count"
   #   >&2 echo "- Count after tagging: \$end_count"
   #   exit 1
   #fi

   mv "umicollapse.bam" "umicollapse-prep.bam"
   samtools view -h "umicollapse-prep.bam" \
   | awk -F'\\t' -v OFS='\\t' '
      /^@/ { print \$0; next }
      !/^@/ {
         delete a; a2l=""; a3l="";
         split(\$1, a, "_"); 
         a2l = length(a[2]);
         a3l = length(a[3]);
         \$1 = a[1] "_" a[2] "_" substr(a[3], a2l + 1, a3l - a2l); 
         print \$0
         next
      }
   ' \
   | samtools view -h -@${task.cpus} -bS -o "umicollapse.bam"

   rm "umicollapse-prep.bam"

   """
}
