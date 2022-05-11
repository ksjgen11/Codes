

#set the directory paths
chippath=/NextGenSeqData/project-data/sejun/chipseq/
webpath=$chippath/web/
httppath=https://bricweb.sund.ku.dk/bric-data/sejun/chipseq/
linkpath=/webdata/bric-data/sejun/chipseq/

# execute BigWig files
makeBigWig.pl $chippath/2/CTCF-ab/ hg38 -normal -force -url $httppath -webdir $linkpath
makeBigWig.pl $chippath/2/IgG/ hg38 -normal -force -url $httppath -webdir $linkpath
makeBigWig.pl $chippath/2/Input/ hg38 -normal -force -url $httppath -webdir $linkpath


