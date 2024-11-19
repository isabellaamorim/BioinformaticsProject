#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <filesystem>

using namespace std;

// Função para converter DNA em RNA
void fromDnaToRna(const string& dna, string& rna) {
    rna.resize(dna.size()); // Redimensiona 'rna' para ter o mesmo tamanho que 'dna'

    #pragma omp parallel for
    for (size_t i = 0; i < dna.size(); ++i) {
        switch (dna[i]) {
            case 'A': rna[i] = 'U'; break;
            case 'T': rna[i] = 'A'; break;
            case 'C': rna[i] = 'G'; break;
            case 'G': rna[i] = 'C'; break;
            default: rna[i] = 'N'; break; // Atribui 'N' a caracteres inválidos
        }
    }
}

// Função para processar arquivos FASTA e chamar a função de conversão
void processFASTA(const string& filename, string& rnaSequence) {
    ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo: " << filename << std::endl;
        return;
    }

    string line, dnaSequence;

    while (getline(inputFile, line)) {
        if (line.empty() || line[0] == '>') continue; // Ignorar linhas de descrição
        dnaSequence += line;
    }

    inputFile.close();

    fromDnaToRna(dnaSequence, rnaSequence);
}

// Função para escrever a sequência de RNA em um arquivo
void writeRNA(const string& filename, const string& rnaSequence) {
    ofstream outputFile(filename);
    if (!outputFile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo: " << filename << std::endl;
        return;
    }

    outputFile << rnaSequence;
    outputFile.close();
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

    // Nº de arquivos por processo
    int numFiles = fastaFiles.size();
    int filesPerProcess = (numFiles + size - 1) / size;

    // Arquivos que cada processo processará
    int start = rank * filesPerProcess;
    int end = min(start + filesPerProcess, numFiles);

    // Processar arquivos FASTA do processo atual
    for (int i = start; i < end; i++) {
        string rnaSequence;
        processFASTA(fastaFiles[i], rnaSequence);

        size_t lastindex = fastaFiles[i].find_last_of(".");

        filesystem::path inputPath(fastaFiles[i]);
        string outputFilename = inputPath.stem().string() + "_rna.txt";

        writeRNA(outputFilename, rnaSequence);

        cout << "Processo " << rank << " converteu o arquivo " << fastaFiles[i] << " para RNA." << endl;
    }

    MPI_Finalize();
    return 0;
}