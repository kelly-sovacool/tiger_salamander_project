#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kellysovacool@uky.edu
SBATCH_NODELIST: $SBATCH_NODELIST

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/all-snps-all-loci_strx.txt -K 7 -L 4948 -N 362.0 -o all-snps-all-loci_k7_run1.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/all-snps-all-loci_strx.txt -K 7 -L 4948 -N 362.0 -o all-snps-all-loci_k7_run2.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/all-snps-all-loci_strx.txt -K 7 -L 4948 -N 362.0 -o all-snps-all-loci_k7_run3.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/all-snps-all-loci_strx.txt -K 7 -L 4948 -N 362.0 -o all-snps-all-loci_k7_run4.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/all-snps-all-loci_strx.txt -K 7 -L 4948 -N 362.0 -o all-snps-all-loci_k7_run5.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/all-snps-all-loci_strx.txt -K 7 -L 4948 -N 362.0 -o all-snps-all-loci_k7_run6.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/all-snps-all-loci_strx.txt -K 7 -L 4948 -N 362.0 -o all-snps-all-loci_k7_run7.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/all-snps-all-loci_strx.txt -K 7 -L 4948 -N 362.0 -o all-snps-all-loci_k7_run8.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/all-snps-all-loci_strx.txt -K 7 -L 4948 -N 362.0 -o all-snps-all-loci_k7_run9.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/all-snps-all-loci_strx.txt -K 7 -L 4948 -N 362.0 -o all-snps-all-loci_k7_run10.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/all-snps-all-loci_strx.txt -K 7 -L 4948 -N 362.0 -o all-snps-all-loci_k7_run11.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/all-snps-all-loci_strx.txt -K 7 -L 4948 -N 362.0 -o all-snps-all-loci_k7_run12.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/all-snps-all-loci_strx.txt -K 7 -L 4948 -N 362.0 -o all-snps-all-loci_k7_run13.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/all-snps-all-loci_strx.txt -K 7 -L 4948 -N 362.0 -o all-snps-all-loci_k7_run14.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/all-snps-all-loci_strx.txt -K 7 -L 4948 -N 362.0 -o all-snps-all-loci_k7_run15.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/all-snps-all-loci_strx.txt -K 7 -L 4948 -N 362.0 -o all-snps-all-loci_k7_run16.out &
sleep 10

wait