
#for i in {34..41}
#do
#    /media/nfs/nfs02/liudy/biosoft/sratoolkit.3.0.2-centos_linux64/bin/fastq-dump --split-3 --gzip -A SRR96670${i} -O ./SRR96670${i}/
#done
#for i in {50..57}
#do
#    /media/nfs/nfs02/liudy/biosoft/sratoolkit.3.0.2-centos_linux64/bin/fastq-dump --split-3 --gzip -A SRR96670${i} -O ./SRR96670${i}/
#done
#for i in {55..62}
#do
#    /media/nfs/nfs02/liudy/biosoft/sratoolkit.3.0.2-centos_linux64/bin/fastq-dump --split-3 --gzip -A SRR117702${i} -O ./SRR117702${i}/
#done
#SRR2339350
#for i in {34..41}
#do
#  /media/nfs/nfs02/wangyi/software/sratoolkit/sratoolkit.3.1.1-centos_linux64/bin/prefetch SRR96670${i}
#done
for i in {405..433}
do
    /media/nfs/nfs02/liudy/biosoft/sratoolkit.3.0.2-centos_linux64/bin/fastq-dump --split-3 --gzip -A SRR5168${i} -O ./SRR5168${i}/
done

