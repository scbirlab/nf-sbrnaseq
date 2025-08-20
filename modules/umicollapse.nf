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
   java -jar -Xmx${Math.round(task.memory.getGiga() * 0.8)}G -Xss1024m \
      ${umicollapse_repo}/umicollapse.jar bam \
      --paired \
      --tag \
      -i ${bamfile[0]} \
      -o ${id}.umicollapse.bam \
   2>&1 > ${id}.umicollapse.log
   """
}
