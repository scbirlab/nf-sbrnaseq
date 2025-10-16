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

process UMIcollapse {

   tag "${id}" 

   label 'big_mem'
   // errorStrategy 'retry'
   // maxRetries 2

   publishDir( 
      "${params.outputs}/umicollapse", 
      mode: 'copy',
   )
   
   input:
   tuple val( id ), path( bamfile ), path( umicollapse_repo )
   val paired

   output:
   tuple val( id ), path( "*.umicollapse.bam" ), emit: main
   path "*.log", emit: logs

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
   ' \
   | samtools sort - -@${task.cpus} -m 2G -t CB \
   | samtools view -bS - -o "cb-prepended.bam"

   mkdir -p shards
   samtools view -H cb-prepended.bam > shards/_hdr.sam

   samtools view -@${task.cpus} cb-prepended.bam \
   | awk -v dir="shards" '
      BEGIN { FS = OFS = "\\t" }
      {
         cb="";
         for(i=12;i<=NF;i++) if(\$i ~ /^CB:Z:/){ split(\$i,a,":"); cb=a[3]; break }
         if(cb == "") next
         pref = substr(cb,1,4); 
         if(pref == "") pref = "NNNN"
         print >> (dir "/" pref ".sam")
      }
   '

   for f in shards/*.sam; do
      if [ "\$f" != "shards/_hdr.sam" ]
      then
         b=\${f%.sam}.bam
         cat shards/_hdr.sam "\$f" \
         | samtools view -@${task.cpus} -b -o "\$b" -
         rm -f "\$f"

         # sort by coordinate and stamp header as coordinate
         samtools sort -@${task.cpus} -m 2G -o "\$b".coord.bam "\$b"
         # ensure @HD SO:coordinate (reheader if needed)
         samtools view -H "\$b".coord.bam \
         | awk '
            BEGIN { done = 0 } 
            /^@HD/ && \$0 ~ /SO:/ { sub(/SO:[^ \\t]+/,"SO:coordinate"); done=1 } 
            { print } 
            END { if(!done) print "@HD\\tVN:1.6\\tSO:coordinate" }
         ' \
         > hdr.coord.sam
         samtools reheader hdr.coord.sam "\$b".coord.bam > "\$b".coord.SOcoord.bam
         samtools index "\$b".coord.SOcoord.bam

         java -jar \
            -Xmx${Math.round(task.memory.getGiga() * 0.8)}G \
            -Xss1024m \
               "${umicollapse_repo}/umicollapse.jar" \
               bam ${paired ? "--paired" : ""} \
               --tag \
               -i "\$b".coord.SOcoord.bam \
               -o "\$b".dedup.bam \
         2>&1 >> "${id}.umicollapse.log"
         rm "\$b".coord.SOcoord.bam hdr.coord.sam "\$b".coord.bam
      fi
   done

   if [ ! -e "${id}.umicollapse.log" ]
   then
      echo "There were no reads in ${bamfile[1]}" > "${id}.umicollapse.log"
   fi

   if ls shards/*.dedup.bam > /dev/null 2>&1
   then
      samtools merge -@${task.cpus} \
         -h cb-prepended.bam \
         -o merged.dedup.bam \
         shards/*.dedup.bam
   else
      samtools view -H cb-prepended.bam -o merged.dedup.bam
   fi

   samtools sort -@${task.cpus} -m 2G \
      -o "${id}.umicollapse.bam" \
      merged.dedup.bam

   samtools index "${id}.umicollapse.bam"
   rm merged.dedup.bam
   
   start_count=\$(samtools view -c "cb-prepended.bam")
   end_count=\$(samtools view -c "${id}.umicollapse.bam")
   if [ "\$start_count" -gt "\$end_count"]
   then
      >&2 echo "Lost some reads during duplicate tagging!"
      >&2 echo "- Initial count: \$start_count"
      >&2 echo "- Count after tagging: \$end_count"
      exit 1
   fi

   mv "${id}.umicollapse.bam" "${id}.umicollapse-prep.bam"
   samtools view -h "${id}.umicollapse-prep.bam" \
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
   | samtools view -h -bS -o "${id}.umicollapse.bam"

   rm "cb-prepended.bam" "${id}.umicollapse-prep.bam"

   """
}
