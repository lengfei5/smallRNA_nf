###################
# count the spkie-ins sequences for small RNA-seq data (initially)
###################
bam="$1"
count=''
umi=''

for seq in Spike_in_X{1..8}
do
    cc=`samtools view ${bam} | grep ${seq} | wc -l`
    uu=`samtools view ${bam} | grep ${seq} | cut -f1|tr '_' '\t'|cut -f2,3|tr '\t' '_'|sort -u |wc -l`
    count="${count} ${cc}"
    umi="${umi} ${uu}"
done
echo spikeIn_{1..8} spikeIn_{1..8}.UMI |tr ' ' '\t'
echo ${count} ${umi} | tr ' ' '\t'
