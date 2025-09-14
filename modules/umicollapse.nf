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

   output:
   tuple val( id ), path( "*.bam" ), emit: main
   path( "*.log" ), emit: logs

   script:
   """
   samtools view -h "${bamfile[1]}" \
   | awk -F'\\t' -v OFS='\\t' '
      /^@/ { print \$0; next }
      !/^@/ {
         delete a
         split(\$1, a, "_"); 
         \$1 = a[1] "_" a[2] "_" a[2] a[3]; 
         print \$0, "CB:Z:" a[2];
         next
      }
   ' \
   | samtools view -h -bS -o "cb-prepended.bam"

   samtools index "cb-prepended.bam"
   java -jar -Xmx${Math.round(task.memory.getGiga() * 0.8)}G -Xss1024m \
      ${umicollapse_repo}/umicollapse.jar bam \
      --paired \
      --tag \
      -i "cb-prepended.bam" \
      -o "${id}.umicollapse.bam" \
   2>&1 > "${id}.umicollapse.log"

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
