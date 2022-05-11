index=~/data/rnaseq/salmon_hg38_index/
trimmed=~/data/rnaseq/trimmed/

salmon quant -i $index -l A --seqBias --gcBias --useVBOpt --validateMappings \
-r $trimmed/LD001_trimmed.fa -o LD001_salmon -p 4
salmon quant -i $index -l A --seqBias --gcBias --useVBOpt --validateMappings \
-r $trimmed/LD002_trimmed.fa -o LD002_salmon -p 4
salmon quant -i $index -l A --seqBias --gcBias --useVBOpt --validateMappings \
-r $trimmed/LD005_trimmed.fa -o LD005_salmon -p 4
salmon quant -i $index -l A --seqBias --gcBias --useVBOpt --validateMappings \
-r $trimmed/LD006_trimmed.fa -o LD006_salmon -p 4
salmon quant -i $index -l A --seqBias --gcBias --useVBOpt --validateMappings \
-r $trimmed/LD007_trimmed.fa -o LD007_salmon -p 4
salmon quant -i $index -l A --seqBias --gcBias --useVBOpt --validateMappings \
-r $trimmed/LD008_trimmed.fa -o LD008_salmon -p 4
salmon quant -i $index -l A --seqBias --gcBias --useVBOpt --validateMappings \
-r $trimmed/LD009_trimmed.fa -o LD009_salmon -p 4
salmon quant -i $index -l A --seqBias --gcBias --useVBOpt --validateMappings \
-r $trimmed/LD010_trimmed.fa -o LD010_salmon -p 4
salmon quant -i $index -l A --seqBias --gcBias --useVBOpt --validateMappings \
-r $trimmed/LD011_trimmed.fa -o LD011_salmon -p 4

