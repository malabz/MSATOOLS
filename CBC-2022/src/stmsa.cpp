#include "StarAlignment/StarAligner.hpp"
#include "Utils/Fasta.hpp"
#include "Utils/Arguments.hpp"
#include <fstream>
#include <string>

int main(int argc, char *argv[])
{
    int zt = 0;
    if (argc==2)
    {
        if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-H") == 0 || strcmp(argv[1], "-help") == 0)
        {
            std::cout << "Usage: softname -ref ref.fasta -reads reads.fastq -out result.txt\n";
            exit(0);
        }
    }
    else if(argc == 7)
    { 
        for (zt = 1; zt < argc; zt+=2)
        {
            if (strcmp(argv[zt], "-ref") == 0)
                arguments::in_file_name = argv[zt + 1];
            else if (strcmp(argv[zt], "-reads") == 0)
                arguments::reads_file_name = argv[zt + 1];
            else if (strcmp(argv[zt], "-out") == 0)
                arguments::out_file_name = argv[zt + 1];
            else
                break;
        }
    }
    else
    {
        std::cout << "Parameter quantity error!\n";
        exit(-1);
    }
    if (zt != argc)
    {
        std::cout << "Parameter error!\n";
        exit(-1);
    }
    std::cout << " ref    name: " << arguments::in_file_name << "\n";
    std::cout << " reads  name: " << arguments::reads_file_name << "\n";
    std::cout << " result name: " << arguments::out_file_name << "\n";
/*
#if defined(_WIN32)
    arguments::in_file_name = "D:/zt/python/code/string/CBC-data/data/";
    arguments::reads_file_name = "D:/zt/python/code/string/CBC-data/data/";
    arguments::out_file_name += "D:/zt/python/code/string/CBC-data/test/";
    
#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    arguments::in_file_name = "/dev/shm/zt/data/CBC/";
    arguments::reads_file_name = "/dev/shm/zt/data/CBC/";
    arguments::out_file_name = "/home/zqzhoutong/CBC/";
#endif  
    arguments::in_file_name += "chr1-3_hg38.fasta";
    arguments::reads_file_name += "simulate_read.sample.unalign.fastq";
    arguments::out_file_name += "align.fasta";
*/
    bool mg_tag = false;

    const auto start_point = std::chrono::high_resolution_clock::now(); //记录总起始时间
    int thresh1 = 15;  //100     全局BWT寻找大于*的同源区间       //两条序列长度，相似性，
    int thresh2 = 1000; //10000    未必对区域分割到*以下                  //根据分割出来的未必对区间，最大长度，来确定
    int thresh3 = 50;    //100  记录大于*的block到.maf结果中,            //用户指定
   
    
    std::vector<std::string> name;
    std::vector<std::size_t> all_Length;
    std::vector<std::size_t> Length;
    std::vector<std::vector<utils::Insertion>> N_gap;
    std::vector<std::vector<unsigned char>> pseudo_sequences;
    if ((arguments::in_file_name.substr(arguments::in_file_name.find_last_of('.') + 1).compare("fasta")) == 0)
    {
        std::ifstream ifs(arguments::in_file_name, std::ios::binary | std::ios::in); //判断输入路径合法否
        if (!ifs)
        {
            std::cout << "cannot access file " << arguments::in_file_name << '\n';
            exit(0);
        }
        utils::read_to_pseudo(pseudo_sequences, ifs, all_Length, Length, name, N_gap);  //读入数据
        ifs.clear();
    }
    else
    {
        std::cout << "The reference sequence format is not .fasta\n";
        exit(-1);
    }

    //std::cout <<name.size()<<" " << name[0] << " " << name[1] << "\n";
    std::cout << Length.size() << " ref_sequences found\n";
    //std::cout << pseudo_sequences.size() << " - ref_sequences size\n";
    
    const auto sy_start = std::chrono::high_resolution_clock::now(); //记录建立索引起始时间
    std::vector<suffix_tree::SuffixTree<nucleic_acid_pseudo::NUMBER> *> bwta(Length.size());
    for (int i = 0; i < Length.size(); i++)
    {
        bwta[i] = new suffix_tree::SuffixTree<nucleic_acid_pseudo::NUMBER>(pseudo_sequences[i].cbegin(), pseudo_sequences[i].cend(), nucleic_acid_pseudo::GAP);
        std::cout << i+1<<" reference sequence BWT is finished.\n";
    }

    std::cout << "memory usage: " << getPeakRSS() << " B" << std::endl;//输出内存耗费
    std::cout << "BWT consumes " << (std::chrono::high_resolution_clock::now() - sy_start) << '\n'; //输出比对耗费时间

    const auto all_start = std::chrono::high_resolution_clock::now(); //记录起始时间
    if ((arguments::reads_file_name.substr(arguments::reads_file_name.find_last_of('.') + 1).compare("fastq")) == 0)
    {
        arguments::fsataq = false;
    }
    else
    {
        std::cout << "The reads sequence format is not .fastq\n";
        exit(-1);
    }
    std::ifstream rfs(arguments::reads_file_name, std::ios::binary | std::ios::in); //判断输入reads路径合法否
    if (!rfs)
    {
        std::cout << "cannot access file " << arguments::reads_file_name << '\n';
        exit(0);
    }
    std::ofstream ofs(arguments::out_file_name, std::ios::binary | std::ios::out); //判断输出路径合法否
    if (!ofs)
    {
        std::cout << "cannot write file " << arguments::out_file_name << '\n';
        exit(0);
    }
    const auto align_start = std::chrono::high_resolution_clock::now(); //记录比对起始时间
    star_alignment::StarAligner::main_pairwise_align(pseudo_sequences, bwta, rfs, ofs, N_gap, name, all_Length, Length, thresh1, thresh2, thresh3);   //****MSA比对过程，insertions记录比对结果
    std::cout << "\n**** process finished****" << std::endl;
    std::cout << "aligning consumes " << (std::chrono::high_resolution_clock::now() - align_start) << '\n'; //输出比对耗费时间
    std::cout << "memory usage: " << getPeakRSS() << " B" << std::endl;//输出内存耗费
    rfs.close();
    ofs.close();
    std::cout << "total consumes" << (std::chrono::high_resolution_clock::now() - all_start) << '\n';//输出总耗费时间
    std::cout << "current pid: " << getpid() << std::endl;
    return 0;
}
