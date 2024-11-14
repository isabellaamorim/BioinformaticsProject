#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

using namespace std;

struct BaseCounts {
    int A = 0, T = 0, C = 0, G = 0;
};

//Contar as bases com OpenMP

void countBases(const std::string& sequence, BaseCounts& counts) {
    int functionA = 0, functionT = 0, functionC = 0, functionG = 0;

    #pragma omp parallel for reduction(+:functionA,functionT,functionC,functionG)
    for (size_t i = 0; i < sequence.size(); ++i) {
        switch (sequence[i]) {
            case 'A': functionA++; break;
            case 'T': functionT++; break;
            case 'C': functionC++; break;
            case 'G': functionG++; break;
        }
    }

    #pragma omp critical
    {
        counts.A += functionA;
        counts.T += functionT;
        counts.C += functionC;
        counts.G += functionG;
    }
}

// Processar arquivos FASTA
void processFASTA(const std::string& filename, BaseCounts& counts) {
    ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo: " << filename << std::endl;
        return;
    }

    string line;

    while (getline(inputFile, line)) {
        if (line.empty() || line[0] == '>') continue; // Ignorar linhas de descrição
        countBases(line, counts);
    }

    inputFile.close();
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    vector<string> fastaFiles = { 
        "chr1.subst.fa", "chr2.subst.fa", "chr3.subst.fa", "chr4.subst.fa", "chr5.subst.fa",
        "chr6.subst.fa", "chr7.subst.fa", "chr8.subst.fa", "chr9.subst.fa", "chr10.subst.fa",
        "chr11.subst.fa", "chr12.subst.fa", "chr13.subst.fa", "chr14.subst.fa", "chr15.subst.fa",
        "chr16.subst.fa", "chr17.subst.fa", "chr18.subst.fa", "chr19.subst.fa", "chr20.subst.fa",
        "chr21.subst.fa", "chr22.subst.fa"
    };

    int numFiles = fastaFiles.size();
    int filesPerProcess = (numFiles + size - 1) / size;

    BaseCounts localCounts;
    BaseCounts globalCounts = {0, 0, 0, 0};

    int start = rank * filesPerProcess;
    int end = min(start + filesPerProcess, numFiles);

    for (int i = start; i < end; i++) {
        processFASTA(fastaFiles[i], localCounts);
    }

    MPI_Reduce(&localCounts.A, &globalCounts.A, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&localCounts.T, &globalCounts.T, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&localCounts.C, &globalCounts.C, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&localCounts.G, &globalCounts.G, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << "Total de bases:" << endl;
        cout << "A: " << globalCounts.A << endl;
        cout << "T: " << globalCounts.T << endl;
        cout << "C: " << globalCounts.C << endl;
        cout << "G: " << globalCounts.G << endl;
    }

    MPI_Finalize();
    return 0;
}