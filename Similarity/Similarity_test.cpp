#include <cstring>
#include <vector>
#include <iostream>

#include "Similarity.h"

void main()
{
    // Test case 1
    unsigned char seq1[] = "AACCATTACTACAGTCGAGCTAATTAGTGCGCAACAAATAGCTTAGTTCC";
    unsigned char seq2[] = "TTCCATTACTACAGTCGAGCTATTAGTGCGCTACAAATAGCATAGTCCC";
    unsigned char seq3[] = "GGAACTAAGCTATTTGTTGCGCACTAATTAGCTCGACTGTAGTAATGGTT";
    size_t len1 = sizeof(seq1) - 1;
    size_t len2 = sizeof(seq2) - 1;
    size_t len3 = sizeof(seq3) - 1;
    int K = 5;
    double sim1 = CaculateSimilarity(K, seq1, seq2, len1, len2);
    std::cout << "Similarity between seq1 and seq2: " << sim1 << std::endl;
    double sim2 = CaculateSimilarity(K, seq1, seq3, len1, len3);
    std::cout << "Similarity between seq1 and seq3: " << sim2 << std::endl;
    if (sim1 > sim2)
        std::cout << "seq1 and seq2 are more similar\n";
    else
        std::cout << "seq1 and seq3 are more similar\n";
}
