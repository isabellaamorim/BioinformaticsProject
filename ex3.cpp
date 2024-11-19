#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

using namespace std;

// Contar AUG utilizando OpenMP
int countAUG(const std::string& sequence) {
    int aug_count = 0;
    #pragma omp parallel for reduction(+:aug_count)
    for (size_t i = 0; i < sequence.size() - 2; ++i) {
        if (sequence[i] == 'A' && sequence[i + 1] == 'U' && sequence[i + 2] == 'G') {
            aug_count++;
        }
    }
    return aug_count;
}

int processFASTA(const std::string& filename) {
    int aug_count = 0;
    ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo: " << filename << std::endl;
        return 0;
    }

    string sequence;
    getline(inputFile, sequence); // Lê a única linha com toda a sequência
    inputFile.close();

    aug_count = countAUG(sequence);
    return aug_count;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    vector<string> fastaFiles = { 
        "chr1.subst_rna.txt", "chr2.subst_rna.txt", "chr3.subst_rna.txt", "chr4.subst_rna.txt", "chr5.subst_rna.txt",
        "chr6.subst_rna.txt", "chr7.subst_rna.txt", "chr8.subst_rna.txt", "chr9.subst_rna.txt", "chr10.subst_rna.txt",
        "chr11.subst_rna.txt", "chr12.subst_rna.txt", "chr13.subst_rna.txt", "chr14.subst_rna.txt", "chr15.subst_rna.txt",
        "chr16.subst_rna.txt", "chr17.subst_rna.txt", "chr18.subst_rna.txt", "chr19.subst_rna.txt", "chr20.subst_rna.txt",
        "chr21.subst_rna.txt", "chr22.subst_rna.txt"
        
    };

    int numFiles = fastaFiles.size();
    int filesPerProcess = (numFiles + size - 1) / size;

    int start = rank * filesPerProcess;
    int end = min(start + filesPerProcess, numFiles);

    int local_aug = 0;
    int global_aug = 0;

    for (int i = start; i < end; i++) {
        local_aug += processFASTA(fastaFiles[i]);
    }

    MPI_Reduce(&local_aug, &global_aug, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "Total de proteínas iniciadas: " << global_aug << endl;
    }

    MPI_Finalize();
    return 0;
}
