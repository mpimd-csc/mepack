#pragma once

#include <tuple>
#include <string>

#ifdef MEXOCT_MATLAB
#include "mexallocator.hpp"
#endif

namespace mexoct {
    namespace Utils {
        // just std::string with allocations mapped to mxMalloc
        // #ifdef MEXOCT_MATLAB
        /*     using string = std::basic_string<char,
               std::char_traits<char>,
               MexAllocator<char>>; */
        // #else
        using string = std::basic_string<char>;
        // #endif

    }

    // internal helpers
    namespace Details{

        // Get index of the first occurrence of a type T in a tuple.
        // Start with Ind = 0. The result is type_index::index.
        template<int Ind, typename T, typename TupleT>
            struct type_index;

        // match found
        template<int Ind, typename T, typename... Ts>
            struct type_index<Ind, T, std::tuple<T,Ts...>>
            {
                constexpr static int index = Ind;
            };

        // recursion through inheritance
        template<int Ind, typename T, typename U, typename... Ts>
            struct type_index<Ind, T, std::tuple<U,Ts...>>
            : type_index<Ind+1,T,std::tuple<Ts...>>
            {};

        // end
        template <int Ind, typename T>
            struct type_index<Ind, T, std::tuple<>>
            {
                constexpr static int index = -1;
            };

        // Check a type_trait on all elements of a parameter pack.
        // source: https://stackoverflow.com/a/29603896
        template<bool...> struct bool_pack;
        template<bool... bs>
            using all_true = std::is_same<bool_pack<bs..., true>, bool_pack<true, bs...>>;


        // simple implementation of std::index_sequence since this is a C++14 feature.
        // https://stackoverflow.com/a/49672613
        template<int... Inds>
            struct index_sequence {};

        template <int N, int ... Next>
            struct index_sequence_helper : public index_sequence_helper<N-1, N-1, Next...>
        { };

        template <int ... Next>
            struct index_sequence_helper<0, Next ... >
            {
                using type = index_sequence<Next ... >;
            };

        template <int N>
            using make_index_sequence = typename index_sequence_helper<N>::type;
    }
}
