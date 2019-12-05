#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kellysovacool@uky.edu
SBATCH_NODELIST: $SBATCH_NODELIST

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 2 -L 92 -N 362.0 -o snp_sub_1_k2_run1.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 2 -L 92 -N 362.0 -o snp_sub_1_k2_run2.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 2 -L 92 -N 362.0 -o snp_sub_1_k2_run3.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 2 -L 92 -N 362.0 -o snp_sub_1_k2_run4.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 2 -L 92 -N 362.0 -o snp_sub_1_k2_run5.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 2 -L 92 -N 362.0 -o snp_sub_1_k2_run6.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 2 -L 92 -N 362.0 -o snp_sub_1_k2_run7.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 2 -L 92 -N 362.0 -o snp_sub_1_k2_run8.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 2 -L 92 -N 362.0 -o snp_sub_1_k2_run9.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 2 -L 92 -N 362.0 -o snp_sub_1_k2_run10.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 2 -L 92 -N 362.0 -o snp_sub_1_k2_run11.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 2 -L 92 -N 362.0 -o snp_sub_1_k2_run12.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 2 -L 92 -N 362.0 -o snp_sub_1_k2_run13.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 2 -L 92 -N 362.0 -o snp_sub_1_k2_run14.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 2 -L 92 -N 362.0 -o snp_sub_1_k2_run15.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 2 -L 92 -N 362.0 -o snp_sub_1_k2_run16.out &
sleep 10

wait