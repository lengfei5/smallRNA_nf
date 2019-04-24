PATH_SINGULARITY_CONTAINER="/groups/cochella/jiwang/scripts/smallRNA_nf/singularity"
if [ `module list -l |fgrep oracle-jdk/1.7.0|wc -l` -eq 1 ]; then
    module unload oracle-jdk/1.7.0
fi;
module load oracle-jdk/1.8.0_72;

bamDir="/groups/cochella/jiwang/Projects/Chiara/NGS_requests/miRNA_R5922_R6016/ngs_raw/*.bam"
outDir="/groups/cochella/jiwang/Projects/Chiara/NGS_requests/miRNA_R5922_R6016"

# run it
~/local/bin/nextflow run -with-singularity $PATH_SINGULARITY_CONTAINER/smallRNA main.nf \
--outdir "$outDir" --reads "$bamDir" -with-timeline smRNA_Chiara.html 

# resume unfinished one
#~/local/bin/nextflow run -resume -with-singularity $PATH_SINGULARITY_CONTAINER/smallRNA main.nf \
#-with-timeline smRNA_spikeIn_trimUMI.html -with-trace -with-report 
