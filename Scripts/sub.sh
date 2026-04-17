#$ -cwd -V

## Runtime in hh:mm:ss
#$ -l h_rt=1:00:00
#$ -l h_vmem=4G
#$ -m be
#$ -M your@email.com
#$ -N your_job_name

## Number of MPI Tasks
#$ -pe ib 16 

mpirun ./main.out >> OUT
