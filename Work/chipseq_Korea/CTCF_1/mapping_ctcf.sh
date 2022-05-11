
filepath=/mnt/d/work/ChIPseq

bowtie2 -x $filepath/references/hg38/GRCh38_noalt_as -1 $filepath/trimmed/CTCF-ab_1_trimmed.fa -2 $filepath/trimmed/CTCF-ab_2_trimmed.fa -S $filepath/CTCF-ab.sam --no-unal --local -p 32
bowtie2 -x $filepath/references/hg38/GRCh38_noalt_as -1 $filepath/trimmed/Ig-G_1_trimmed.fa -2 $filepath/trimmed/Ig-G_2_trimmed.fa -S $filepath/Ig-G.sam --no-unal --local -p 32
bowtie2 -x $filepath/references/hg38/GRCh38_noalt_as -1 $filepath/trimmed/Input_1_trimmed.fa -2 $filepath/trimmed/Input_2_trimmed.fa -S $filepath/Input.sam --no-unal --local -p 32

