# set the directory paths
chippath=/NextGenSeqData/project-data/sejun/chipseq/
webpath=$chippath/web/
httppath=https://bricweb.sund.ku.dk/bric-data/sejun/chipseq/
linkpath=/webdata/bric-data/sejun/chipseq/CTCF_IgG

# execute multiwighub to create UCSC hub
makeMultiWigHub.pl CTCF_IgG hg38 -d $chippath/2/CTCF-ab/ $chippath/2/Ig-G/ $chippath/2/Input/ -normal -webdir $webpath -url $httppath -force

# make symbolic link in webdata folder
ln -s $webpath/CTCF_IgG/hg38 $linkpath/hg38
ln -s $webpath/CTCF_IgG/hub.txt $linkpath/hub.txt
ln -s $webpath/CTCF_IgG/genomes.txt $linkpath/genomes.txt
