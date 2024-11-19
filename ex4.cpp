#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <mpi.h>
#include <omp.h>

using namespace std;

// Mapeamento de códons para aminoácidos
map<string, char> codonToAmino = {
    {"CCA", 'P'}, {"CCG", 'P'}, {"CCU", 'P'}, {"CCC", 'P'},
    {"UCU", 'S'}, {"UCA", 'S'}, {"UCG", 'S'}, {"UCC", 'S'},
    {"CAG", 'Q'}, {"CAA", 'Q'},
    {"ACA", 'T'}, {"ACC", 'T'}, {"ACU", 'T'}, {"ACG", 'T'},
    {"AUG", 'M'}, 
    {"UAA", 'X'}, {"UAG", 'X'}, {"UGA", 'X'}, // Códons de parada
    {"UGC", 'C'}, {"UGU", 'C'},
    {"GUG", 'V'}, {"GUA", 'V'}, {"GUC", 'V'}, {"GUU", 'V'}
};

// Tradução de RNA para proteínas
string getProteins(const string& rna) {
    string protein;
    for (size_t i = 0; i < rna.size(); i += 3) {
        string codon = rna.substr(i, 3);
        if (codonToAmino.find(codon) != codonToAmino.end()) {
            char amino = codonToAmino[codon];
            protein += amino;

            // Códon de parada adiciona uma nova linha
            if (amino == 'X') {
                protein += '\n';
            }
        } 
    }
    return protein;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    cout << "[INFO] Processo " << rank << " iniciado com " << size << " processos no total." << endl;

    vector<string> rnaFiles = { 
        "chr1.subst_rna.txt", "chr2.subst_rna.txt", "chr3.subst_rna.txt", "chr4.subst_rna.txt", "chr5.subst_rna.txt",
        "chr6.subst_rna.txt", "chr7.subst_rna.txt", "chr8.subst_rna.txt", "chr9.subst_rna.txt", "chr10.subst_rna.txt",
        "chr11.subst_rna.txt", "chr12.subst_rna.txt", "chr13.subst_rna.txt", "chr14.subst_rna.txt", "chr15.subst_rna.txt",
        "chr16.subst_rna.txt", "chr17.subst_rna.txt", "chr18.subst_rna.txt", "chr19.subst_rna.txt", "chr20.subst_rna.txt",
        "chr21.subst_rna.txt", "chr22.subst_rna.txt"
    };

    // Divisão dos arquivos entre os processos
    int numFiles = rnaFiles.size();
    int filesPerProcess = (numFiles + size - 1) / size;

    int start = rank * filesPerProcess;
    int end = min(start + filesPerProcess, numFiles);

    cout << "[INFO] Processo " << rank << " processará arquivos de " << start << " até " << end - 1 << endl;

    // Processo atual
    vector<string> localResults;

    // Processamento dos arquivos
    for (int i = start; i < end; i++) {
        ifstream inputFile(rnaFiles[i]);
        if (!inputFile.is_open()) {
            cerr << "[ERROR] Erro ao abrir o arquivo: " << rnaFiles[i] << endl;
            continue;
        }

        string rna;
        getline(inputFile, rna);
        inputFile.close();
        cout << "[INFO] Arquivo " << rnaFiles[i] << " lido com sucesso pelo processo " << rank << endl;

        string translated;
        #pragma omp parallel
        {
            #pragma omp single
            translated = getProteins(rna);
        }

        localResults.push_back(translated);
        cout << "[INFO] RNA traduzido do arquivo " << rnaFiles[i] << " pelo processo " << rank << endl;
    }

    // Agregação dos resultados
    ofstream outFile("output_" + to_string(rank) + ".txt");
    if (!outFile.is_open()) {
        cerr << "[ERROR] Não foi possível criar o arquivo de saída para o processo " << rank << endl;
    } else {
        for (const auto& result : localResults) {
            outFile << result << endl;
        }
        outFile.close();
    }

    MPI_Finalize();
    cout << "[INFO] Processo " << rank << " finalizado." << endl;

    return 0;
}
