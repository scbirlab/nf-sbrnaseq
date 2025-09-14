// Build whitelist from provided barcode files
process build_whitelist {

   tag "${id}:${reverse ? "rev ${reverse}" : ""}"

   publishDir( 
      "${params.outputs}/barcodes", 
      mode: 'copy',
      pattern: "whitelist.txt.gz",
      saveAs: { "${id}-${it}" },
   )

   input:
   tuple val( id ), path( bcs )
   val reverse

   output:
   tuple val( id ), path( "whitelist.txt.gz" )

   script:
   """
   set -euox pipefail
   BARCODES=(${bcs})
   for i in "\${!BARCODES[@]}"
   do
      i_p1=\$(("\$i"+1))
      f="\${BARCODES[\$i]}"

      tr -d \$'\\r' \
      < "\$f" \
      > "\$f".clean \
      && mv "\$f".clean "\$f"

      cp "\$f" "\$f".inp

      for j in ${reverse}
      do
         if [ "\$i_p1" -eq "\$j" ]
         then
            paste -d, \
             <(
               tail -n+2 "\$f" \
               | cut -f1 -d,
             ) \
             <(
               tail -n+2 "\$f" \
               | cut -f2 -d, \
               | tr ATCGNatcgn TAGCNtagcn \
               | rev
             ) \
            > "\$f".rev

            head -n1 "\$f" \
            | cat - "\$f".rev \
            > "\$f".inp
         fi
      done

      tail -n+2 "\$f".inp \
      | awk -F, -v OFS=, '{ print "__JOIN__", \$2 }' \
      > "\$f.temp"
   done

   join -t, ${bcs[0]}.temp ${bcs[1]}.temp \
   | join -t, - ${bcs[2]}.temp \
   | awk -F, '{ print \$2\$3\$4 }' \
   > whitelist.txt

   pigz -v -p ${task.cpus} whitelist.txt

   """
}

// Build whitelist from provided barcode files
process add_errors_to_whitelist {

   tag "${id}"
   label "big_time"

   publishDir( 
      "${params.outputs}/barcodes", 
      mode: 'copy',
      pattern: "whitelist-err.txt.gz",
      saveAs: { "${id}-${it}" },
   )

   input:
   tuple val( id ), path( whitelist )

   output:
   tuple val( id ), path( "whitelist-err.txt.gz" )

   script:
   """
   #!/usr/bin/env python

   from itertools import product
   import gzip

   ALPHABET = "ATCGN"
   INPUT_FILE = "${whitelist}"
   OUTPUT_FILE = "whitelist-err.txt.gz"

   with gzip.open(OUTPUT_FILE, 'wt') as outfile, gzip.open(INPUT_FILE, 'rt') as infile:
      for line in infile:
         bc_length = len(line)
         break
      infile.seek(0)
      bc_too_long = (bc_length * (len(ALPHABET) - 1)) > 120
      if not bc_too_long:
         for line in infile:
            line = line.strip()
            alternatives = (
               f"{line[:i]}{letter}{line[(i+1):]}" 
               for (i, char), letter in product(
                  enumerate(line), 
                  ALPHABET,
               ) if char != letter
            )
            alternatives = ','.join(sorted(set(alternatives)))
            print(f"{line}\\t{alternatives}", file=outfile)
      else:
         for line in infile:
            print(line.strip(), file=outfile)

   """
}


process chunk_whitelist {

   tag "${id}"
   label "big_cpu"

   // publishDir( 
   //    "${params.outputs}/barcodes", 
   //    mode: 'copy',
   //    pattern: "whitelist.txt",
   //    saveAs: { "${id}-${it}" },
   // )

   input:
   tuple val( id ), path( whitelist )
   val chunk_size

   output:
   tuple val( id ), path( "chunk_*.gz" )

   script:
   """
   zcat "${whitelist}" | split -l "${chunk_size}" - chunk_
   for f in chunk_*
   do
      pigz -v -p ${task.cpus} "\$f"
   done

   """
}

