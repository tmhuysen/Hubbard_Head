//
// Created by Wulfix on 17/11/2017.
//

#ifndef HUBBARDCLEAN_BITSET_H
#define HUBBARDCLEAN_BITSET_H

#include <boost/dynamic_bitset.hpp>
#include <boost/math/special_functions.hpp>

/** Give the next bitset permutation in a lexicographical sense.
 *
 *      examples:
 *          011 -> 101
 *          101 -> 110
 *
 * Taken from (http://graphics.stanford.edu/~seander/bithacks.html#NextBitPermutation)
 * Credits: Dario Sneidermanis of Argentina, who provided this on November 28, 2009.
 */
boost::dynamic_bitset<> next_bitset_permutation(const boost::dynamic_bitset<>& v);



/** Generate all bitset permutations in a lexicographical sense, with the lowest being given as an argument.
 *
 *      example:
 *          011 -> std::vector {011, 101, 110}
 */
std::vector<boost::dynamic_bitset<>> all_bitset_permutations(const boost::dynamic_bitset<>& v);


/** Generate the smallest (in lexicographical sense) bitset, given a total length and number of ones (nones).
 *
 *      example:
 *          length=3, nones=2 -> 011
 */
boost::dynamic_bitset<> smallest_bitset(size_t length, size_t nones);


/** Checks if the given bitset is the lexicographical smallest one.
 */
bool is_smallest_bitset(const boost::dynamic_bitset<>& v);


/** Concatenate two sets of bitsets, producing every possible combination between both. set1 is major.
 *
 *      example:
 *          concatenate_bitset_sets {000}, {011, 101, 110} -> {000011, 000101, 000110}
 *
 *
 * @param set1
 * @param set2
 * @return
 */
std::vector<boost::dynamic_bitset<>> concatenate_bitset_sets(const std::vector<boost::dynamic_bitset<>>& set1, const std::vector<boost::dynamic_bitset<>>& set2);


/** Count the number of differences between two bitsets.
 *
 *      example:
 *          count_differences(001, 100) -> 2
 *          count_differences(001, 011) -> 1
 */
unsigned long count_differences(const boost::dynamic_bitset<>& u, const boost::dynamic_bitset<>& v);


/** Slice a given bitset from least significant to most significant bit (i.e. from right to left) This is NOT as they are printed!
 *
 *     example:
 *          u: 0100011101
 *          bitslice(u, 2, 5) is NOT equal to 0001
 *          bitslice(u, 2, 5) => 0100[0111]01 -> 0111
 */
boost::dynamic_bitset<> bitslice(const boost::dynamic_bitset<>& u, unsigned long i, unsigned long j);

/** Count the number of double occupancies between two bitsets
 *
 *      example:
 *          count_double_occupancies(00110, 10101) -> 1
 *          count_double_occupancies(1110, 1100) -> 2
 */
unsigned long count_double_occupancies(const boost::dynamic_bitset<>& u, const boost::dynamic_bitset<>& v);


/**
 *  Applies the annihilation operator (without sign evaluation)
 *  returns bool according to succes.
 */
bool annihilation( boost::dynamic_bitset<>& u, unsigned long i);
/**
 *  Applies the creation operator (without sign evaluation)
 *  returns bool according to succes.
 */
bool creation( boost::dynamic_bitset<>& u, unsigned long i);

/**
 * Evaluates the sign change between two (spin spererated) base vectors
 * (only for vectors changed by one creation and one annihilation operators)
 */

int phaseCheck( const boost::dynamic_bitset<>& u, const boost::dynamic_bitset<>& v);

#endif //HUBBARDCLEAN_BITSET_H
