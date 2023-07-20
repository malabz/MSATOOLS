#pragma once
//.fasta��ȡ��д��
#include <string>
#include <vector>
#include <iostream>

namespace utils
{

    class Fasta
    {
    private:
        void _read(std::istream &is);//��ȡ�ļ�

    public:
        static constexpr unsigned max_line_length = 80; //һ���

        std::vector<std::string> sequences;
        std::vector<std::string> identifications;

        explicit Fasta(std::istream &is);

        void write_to(std::ostream &os, bool with_idification = true) const;

        static void cut_and_write(std::ostream &os, const std::string &sequence);

        template<typename InputIterator>
        static void write_to(std::ostream &os, InputIterator sequence_first, InputIterator sequence_last)
        {
            if (sequence_first == sequence_last) return;

            using difference_type = decltype(std::distance(sequence_first, sequence_last));
            const difference_type len = std::distance(sequence_first, sequence_last);

            for (difference_type i = 0; i != len; ++sequence_first, ++i)
            {
                os << *sequence_first;
                if (i != len - 1) os << '\n';
            }
        }

        template<typename InputIterator1, typename InputIterator2>
        static void write_to(std::ostream &os, InputIterator1 sequence_first, InputIterator1 sequence_last, //д���ļ�
                InputIterator2 identification_first)
        {
            if (sequence_first == sequence_last) return;

            using difference_type = decltype(std::distance(sequence_first, sequence_last));
            const difference_type len = std::distance(sequence_first, sequence_last);

            for (difference_type i = 0; i != len; ++sequence_first, ++identification_first, ++i)
            {
                os << '>' << *identification_first << '\n';
                cut_and_write(os, *sequence_first);
                if (i != len - 1) os << '\n';
            }
        }

    };

}
