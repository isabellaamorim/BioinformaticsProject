#!/bin/bash
#SBATCH --job-name=ex4               # Nome do job
#SBATCH --output=ex4_output_%j.txt   # Nome do arquivo de saída (%j insere o job ID)
#SBATCH --error=ex4_error_%j.txt     # Nome do arquivo de erro (%j insere o job ID)
#SBATCH --ntasks=8                   # Número total de tarefas MPI (processos)
#SBATCH --nodes=1                    # Número de nós (1 nó para 4 tarefas)
#SBATCH --mem=1024                   # Memória total alocada por nó (1024 MB)
#SBATCH --time=00:05:00              # Tempo limite de execução (hh:mm:ss)
#SBATCH --partition=espec            # Nome da partição/filas 

# Executar o programa com MPI
mpirun -np 8 ./ex4