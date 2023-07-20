#pragma once
#include "../SuffixTree/SuffixTree.hpp"  //利用后缀树
#include "../Utils/Utils.hpp"

#include <vector>
#include <array>
#include <string>
#include<algorithm>


namespace star_alignment //星比对命名空间
{

    
    class StarAligner//星比对类
    {
    private:
        using triple = std::array<size_t, 3>;
        using quadra = std::array<size_t, 4>;
        using sequence_type = std::vector<unsigned char>;

    public:
        struct Substrings
        {
            size_t id;
            bool sign;
            std::vector<std::vector<triple>> Subs;
            size_t add_len;
        };
        struct Sub_align
        {
            size_t id;
            bool sign;
            size_t a_start, b_start, a_end, b_end; //闭区间
            size_t nw_start, nw_end;
        };
        static bool cmp_Substrings(Substrings x, Substrings y);
        static bool cmp_Sub_align(Sub_align x, Sub_align y);
        static void get_common_substrings_vector(std::vector<Substrings>& C_Strings, std::vector<triple>& substrings, int id, bool sign, size_t read_len);
        static void main_pairwise_align(
            std::vector<std::vector<unsigned char>>& _a_sequence,
            std::vector<suffix_tree::SuffixTree<nucleic_acid_pseudo::NUMBER>*>& _bwta,
            std::ifstream& _is,
            std::ofstream& _os,
            std::vector<std::vector<utils::Insertion>>& _N_gap,
            std::vector<std::string>& _name,
            std::vector<std::size_t>& _all_Length,
            std::vector<std::size_t>& _Length,
            int threshold1, int threshold2, int threshold3);

        //依据动态规划选择最长不重叠同源区段
        //传入全部的同源区段，[[A.index，B.index，length]...] B.index可能相同即[[[A.index0,A.index1...],B.index,length],...]
        //返回选择好的同源区段[[A.index，B.index，length]...] B.index各不相同即[[A.index,B.index,length],...]
        static std::vector<triple> _optimal_path(const std::vector<triple> &common_substrings);//最优路径
        //分两步新改,第一步，从传入全部的同源区段，[[A.index，B.index，length]...]，每个B.index选出一个A.index，相对距离最近
        static std::vector<triple> _MultiReg(const std::vector<triple>& common_substrings);
        static std::vector<int> _trace_back(const std::vector<triple>& common_substrings, int* p);
    private:
        StarAligner(
            std::vector<std::vector<unsigned char>>& _a_sequence,
            std::vector<suffix_tree::SuffixTree<nucleic_acid_pseudo::NUMBER>*>& _bwta,
            std::ifstream& _is,
            std::ofstream& _os,
            std::vector<std::vector<utils::Insertion>>& _N_gap,
            std::vector<std::string>& _name,
            std::vector<std::size_t>& _all_Length,
            std::vector<std::size_t>& _Length
        );

        void _pairwise_align(int threshold1, int threshold2, int threshold3);
        void NW2(std::vector<Sub_align>& align, Sub_align& align_tmp, utils::reads& now_read);

        std::vector<std::vector<unsigned char>> &a_sequence;
        std::vector<std::string>& name;
        std::vector<std::size_t>& Length;
        std::vector<std::size_t>& all_Length;
        std::vector<std::vector<utils::Insertion>>& N_gap;
        std::vector<suffix_tree::SuffixTree<nucleic_acid_pseudo::NUMBER>*>& bwta;
        std::ifstream& is;
        std::ofstream& os;
    };

}
