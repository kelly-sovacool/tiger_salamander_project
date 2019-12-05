#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kellysovacool@uky.edu
SBATCH_NODELIST: $SBATCH_NODELIST

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 5 -L 92 -N 362.0 -o snp_sub_1_k5_run1.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 5 -L 92 -N 362.0 -o snp_sub_1_k5_run2.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 5 -L 92 -N 362.0 -o snp_sub_1_k5_run3.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 5 -L 92 -N 362.0 -o snp_sub_1_k5_run4.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 5 -L 92 -N 362.0 -o snp_sub_1_k5_run5.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 5 -L 92 -N 362.0 -o snp_sub_1_k5_run6.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 5 -L 92 -N 362.0 -o snp_sub_1_k5_run7.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 5 -L 92 -N 362.0 -o snp_sub_1_k5_run8.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 5 -L 92 -N 362.0 -o snp_sub_1_k5_run9.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 5 -L 92 -N 362.0 -o snp_sub_1_k5_run10.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 5 -L 92 -N 362.0 -o snp_sub_1_k5_run11.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 5 -L 92 -N 362.0 -o snp_sub_1_k5_run12.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 5 -L 92 -N 362.0 -o snp_sub_1_k5_run13.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 5 -L 92 -N 362.0 -o snp_sub_1_k5_run14.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 5 -L 92 -N 362.0 -o snp_sub_1_k5_run15.out &
sleep 10

/home/klso224/tools/structure-2.3.4/structure -i snp_subsamples/snp_sub_1_strx.txt -K 5 -L 92 -N 362.0 -o snp_sub_1_k5_run16.out &
sleep 10

wait