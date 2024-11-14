#include <iostream>
#include <fstream>
#include <string>
#include <cctype>
#include <vector>
#include <omp.h>

void convertToUpperCase(std::string& sequence) {
    for (char& base : sequence) {
        base = std::toupper(base);
    }
}

void processFASTAFile(const std::string& inputFilename, const std::string& outputFilename) {
    std::ifstream inputFile(inputFilename);
    std::ofstream outputFile(outputFilename);

    if (!inputFile.is_open() || !outputFile.is_open()) {
        std::cerr << "Erro ao abrir o arquivo: " << inputFilename << " ou " << outputFilename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        if (line.empty()) continue;

        // Se a linha começar com '>', é uma linha de descrição e não deve ser alterada
        if (line[0] == '>') {
            outputFile << line << std::endl;
        } else {
            convertToUpperCase(line);
            outputFile << line << std::endl;
        }
    }

    inputFile.close();
    outputFile.close();
}

int main() {
    std::vector<std::string> fastaFiles = { 
        "chr1.subst.fa", "chr2.subst.fa", "chr3.subst.fa", "chr4.subst.fa", "chr5.subst.fa",
        "chr6.subst.fa", "chr7.subst.fa", "chr8.subst.fa", "chr9.subst.fa", "chr10.subst.fa",
        "chr11.subst.fa", "chr12.subst.fa", "chr13.subst.fa", "chr14.subst.fa", "chr15.subst.fa",
        "chr16.subst.fa", "chr17.subst.fa", "chr18.subst.fa", "chr19.subst.fa", "chr20.subst.fa",
        "chr21.subst.fa", "chr22.subst.fa"
    };

    // Paraleliza o loop para processar cada arquivo em uma thread separada
    #pragma omp parallel for

    for (size_t i = 0; i < fastaFiles.size(); ++i) {
        const auto& filename = fastaFiles[i];
        std::string outputFilename = "upper_" + filename;
        processFASTAFile(filename, outputFilename);
        
        #pragma omp critical
        {
            std::cout << "Processado: " << filename << " -> " << outputFilename << std::endl;
        }
    }

    return 0;
}
