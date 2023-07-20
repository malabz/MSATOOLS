#pragma once
#include "../SuffixTree/SuffixTree.hpp"  //���ú�׺��
#include "../Utils/Utils.hpp"

#include <vector>
#include <array>
#include <string>
#include<algorithm>


namespace star_alignment //�Ǳȶ������ռ�
{

    
    class StarAligner//�Ǳȶ���
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
            size_t a_start, b_start, a_end, b_end; //������
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

        //���ݶ�̬�滮ѡ������ص�ͬԴ����
        //����ȫ����ͬԴ���Σ�[[A.index��B.index��length]...] B.index������ͬ��[[[A.index0,A.index1...],B.index,length],...]
        //����ѡ��õ�ͬԴ����[[A.index��B.index��length]...] B.index������ͬ��[[A.index,B.index,length],...]
        static std::vector<triple> _optimal_path(const std::vector<triple> &common_substrings);//����·��
        //�������¸�,��һ�����Ӵ���ȫ����ͬԴ���Σ�[[A.index��B.index��length]...]��ÿ��B.indexѡ��һ��A.index����Ծ������
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
