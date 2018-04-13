#include "boost/multi_array.hpp"
#include <boost/math/special_functions/binomial.hpp> //  boost::math::binomial_coefficient<int>(10, 2);
#include <math.h>
#include <vector>
#include "tandem.h"

using namespace std;
Tandem::Tandem(const int num_orbitals, const int num_particles, const vector<double> wave_function)
  : num_particles(num_particles), num_orbitals(num_orbitals),
    wave_function(wave_function), two_rdm_all_zero(false),
    two_rdm(boost::extents[num_orbitals][num_orbitals][num_orbitals]
	    [num_orbitals]) {
    // Initialize resulting 2DRM
    for (int i = 0; i < num_orbitals; ++i) {
      for (int j = 0; j < num_orbitals; ++j) {
        for (int k = 0; k < num_orbitals; ++k) {
          for (int l = 0; l < num_orbitals; ++l) {
            two_rdm[i][j][k][l] = 0;
          }
        }
      }
    }
  };

vector<double> Tandem::run() {
    // If the method has been run before we need to make sure we are
    // starting from a clean slate
    if (not two_rdm_all_zero) {
      for (int i = 0; i < num_orbitals; i++) {
        for (int j = 0; j < num_orbitals; j++) {
          for (int k = 0; k < num_orbitals; k++) {
            for (int l = 0; l < num_orbitals; l++) {
              two_rdm[i][j][k][l] = 0;
            }
          }
        }
      }
    }
    // Allocate our bool arrays. I did not dare use std::vector<bool>
    // as some warnings where made about it by random people on the
    // internet. Maybe there is more problems with this approache but
    // I am trying to take care to not create memmory leaks. And the
    // length is allways the same during one run and there are plenty
    // of tests to make sure it all works. And it have been written
    // and debugged with bitset so I hope this will work out well.
    // So please do not judge, or if you can not stop yourself,
    // do a pull request showing how it should be done!
    
    // These are the bool arrays we will need:
    // ---------------------------------------
    // First two to permute
    bool* state1 = new bool[num_orbitals];
    bool* state2 = new bool[num_orbitals];
    // Then two with nice index to creation operator index convention.
    bool* bitstate1 = new bool[num_orbitals];
    bool* bitstate2 = new bool[num_orbitals];
    // One to know the difference between them
    bool* diff = new bool[num_orbitals];
    // Then two to know where to create and where to anihilate.
    bool* b1 = new bool[num_orbitals];
    bool* b2 = new bool[num_orbitals];
    // Then one to know what is in common
    bool* common = new bool[num_orbitals];
    // That is all that need to be dynammically allocated.
    // Create arrays of bools with num_particles ones
    for (int i = 0; i < num_orbitals; ++i) {
      if (i < num_particles) {
        state1[i] = 1;
        state2[i] = 1;
      } else {
        state1[i] = 0;
        state2[i] = 0;
      }
    }
    // Will switch positions of the 1s so that they are last in the
    // array. Lowest in lexiographical comparison
    sort(state1, state1 + num_orbitals);
    sort(state2, state2 + num_orbitals);
    // We need to keep track on what position in the wave_function that the
    // state corresponds to
    int num_state1 = 0;
    int num_state2 = 0;
    // Iterate through all combinations of 2 permutations of
    // num_particles 1s in a list of length num_orbitals.

    do {
      copy(state2, state2 + num_orbitals, state1); // Start over the state1
                                                   // iteration from state2
      num_state1 = num_state2;
      if (abs(wave_function[num_state2]) > 1e-15) {
        do {
          if (abs(wave_function[num_state1]) > 1e-15) {
            // Easier to make bitwise operations as bitset but bitset
            // does not support iterators so cannot be used with
            // next_permutation The ordering of bits is also changed
            // so instead of counting from the left when translating
            // to creation and anihilation operators we calculate from
            // the right giving the translation as the usual index.
	    
            for (int i = 0; i < num_orbitals; ++i) {
              if (state1[i]) {
                bitstate1[num_orbitals - 1 - i] = 1;
              }
	      else {
		bitstate1[num_orbitals - 1 - i] = 0;
	      }
              if (state2[i]) {
                bitstate2[num_orbitals - 1 - i] = 1;
              }
	      else {
		bitstate2[num_orbitals - 1 - i] = 0;
	      }
            }
            // now that we have nice data format lets get to buissness.
            // Xor to get positions that differ
            elementwise_xor(bitstate1, bitstate2, diff, num_orbitals);
            // Now lets handle this in different cases depending on
            // the number of elements that differ in the two states.
            switch (count(diff, num_orbitals)) {
            case 4: {
              // b1 (b2) - elements that need to be anihilated in bitstate1
              // (bitstate2)
              // for the states not to differ.
              elementwise_and(diff, bitstate1, b1, num_orbitals);
	      elementwise_and(diff, bitstate2, b2, num_orbitals);
              // Whould cause ansty undefined errors if this was not true.
              assert(count(b1, num_orbitals) == 2);
              assert(count(b2, num_orbitals) == 2);
              // The positions where the states differ give whitch index the
              // creation and
              // anhihilation operators must have. (2RDM)_ijkl =
              // <bitstate1|a^d_ia^d_ja_ka_l|bitstate2>
              // So, i, j note the states that need to be
              // anihilitaed in bitstate 1
              // and k,l note the states that need to anhiliated in bitstate2.
              int index[4];
              int counter1 = 0, counter2 = 2;
              for (int i = 0; i < num_orbitals; ++i) {
                if (b1[i]) {
                  index[counter1] = i;
                  counter1++;
                }
                if (b2[i]) {
                  index[counter2] = i;
                  counter2++;
                }
              }
              order_and_save(index, bitstate2, num_state1, num_state2);
              break;
            }
            case (2): {
              // Now we have some more leyway to assign the indexes
              // b1 (b2) - elements that need to be anihilated in bitstate1
              // (bitstate2)
              // for the states not to differ.
              elementwise_and(diff, bitstate1, b1, num_orbitals);
	      elementwise_and(diff, bitstate2, b2, num_orbitals);
              // Whould cause ansty undefined errors if this was not true.
	      assert(count(b1, num_orbitals) == 1);
	      assert(count(b2, num_orbitals) == 1);
              int index[4];
              // Find the position of the 1 and save it to i, and l respectivly
              // counting the position from the end of the bitset
              for (int i = 0; i < num_orbitals; ++i) {
                if (b1[i]) {
                  index[0] = i;
                }
                if (b2[i]) {
                  index[3] = i;
                }
              }
              // Now to the j, k pair. we must have j == k == index of true
              // element in bitstate2 that is not eqal to l.
	      elementwise_and(bitstate1, bitstate2, common, num_orbitals);
              for (int i = 0; i < num_orbitals; ++i) {
                if (common[i]) {
                  index[1] = i;
                  index[2] = i;
                  int temp[4];
                  copy(index, index + 4, temp);
                  order_and_save(temp, bitstate2, num_state1, num_state2);
                }
              }
              break;
            }
            case 0: {
              // the bitstates are the same.  Here all bets are off, total
              // freedom, whoo!! Well not really but a bit more than the
              // other cases. We ensure i == k and j == l, i != j, and
              // bitstate[i] == bitstate[j] == 1.  What we do is that we
              // write down the filled orbitals and combine them in all
              // uniqe combinations
              int non_zero[num_particles];
              int counter = 0;
              for (int i = 0; (i < num_orbitals) && (counter < num_particles);
                   i++) {
                if (bitstate2[i]) {
                  non_zero[counter] = i;
                  counter++;
                }
              }
              sort(non_zero, non_zero + num_particles);
              assert(counter == num_particles);
              for (int i = 0; i < num_particles; i++) {
                for (int j = i + 1; j < num_particles; j++) {
                  int index[4] = {non_zero[i], non_zero[j], non_zero[i],
                                  non_zero[j]};
                  order_and_save(index, bitstate2, num_state1, num_state2);
                }
              }
              break;
            }
            }
            // NOTE: There is no case for count > 4 beacuse these will
            // not give any contributions.
          }
          num_state1++;
        } while (std::next_permutation(state1, state1 + num_orbitals));
      }
      num_state1 = 0;
      num_state2++;
    } while (std::next_permutation(state2, state2 + num_orbitals));
    // Now we have saved everything to two_rdm, that is however a very
    // weird form as not all elements in two_rdm hold any interesting
    // information.  Only the elements with orderd indices hold
    // information therefore we return the thing in the handy .red
    // format as a list.  We use std::vector here as there is a small
    // chance the looping is wrong as it is complicated...
    int basis_length = boost::math::binomial_coefficient<double>(
        num_orbitals, 2); // Two particle basis
    vector<double> res(
        int(boost::math::binomial_coefficient<double>(basis_length, 2) +
            basis_length));
    int counter = 0;
    vector<vector<int>> list_ij(basis_length, vector<int>(2));
    for (int i = 0; i < num_orbitals; i++) {
      for (int j = i + 1; j < num_orbitals; j++) {
        list_ij.at(counter).at(0) = i;
        list_ij.at(counter).at(1) = j;
        counter++;
      }
    }
    counter = 0;
    for (int n = 0; n < list_ij.size(); n++) {
      for (int m = n; m < list_ij.size(); m++) {
        int i = list_ij.at(n).at(0);
        int j = list_ij.at(n).at(1);
        int k = list_ij.at(m).at(0);
        int l = list_ij.at(m).at(1);
        res.at(counter) = two_rdm[i][j][k][l];
        counter++;
      }
    }
    assert(counter == res.size());
    two_rdm_all_zero = false;
    // For every new there must be a delete.  Only using new in one
    // place in the begining of this method and all those things
    // should be deleted here
    delete [] state1;
    delete [] state2;
    delete [] bitstate1;
    delete [] bitstate2;
    delete [] diff;
    delete [] b1;
    delete [] b2;
    delete [] common;
    return res; // Return std::vector<double> in .red matrix order.
  }

void elementwise_xor(bool* x, bool* y, bool* res, int size){
  for (int i = 0; i < size; i++) {
    res[i] = x[i] != y[i];
  }
}

void elementwise_and(bool* x, bool* y, bool* res, int size){
  for (int i = 0; i < size; i++) {
    res[i] = x[i] && y[i];
  }
}

int count(bool* x, int size){
  int c = 0;
  for (int i = 0; i < size; i++){
    if (x[i]){
      c++;
    }
  }
  return c;
}

bool Tandem::sign(int index[4], bool* state) {
    // All this work for a silly sign. But you know what they say: the
    // devil is in the details.
    // check that annihilators have something to anihiliate
    int sign_index[4];
    copy(index, index + 4, sign_index);
    for (int i = 2; i < 4; i++) {
      if (not state[sign_index[i]]) {
        swap(sign_index[0], sign_index[2]);
        swap(sign_index[1], sign_index[3]);
      }
    }

    // Start with index[3], count the number of 1s
    // before position index in state. Actually we just keep track of
    // if it is an even or odd number with a bool.

    // For index[3]/[2], state[index[3/2]] == 1 change that to 0
    // For index[2]/[0], state[index[2/0]] == 0 change that to 1
    bool neg_sign = false; // sign starts as positive.
    for (int i = 3; i >= 0; i = i - 1) {
      state[sign_index[i]] = ! state[sign_index[i]];
      for (int j = sign_index[i] + 1; j < num_orbitals; j++) {
        if (state[j]) {
          neg_sign = ! neg_sign;
        }
      }
    }
    return neg_sign;
  }
void Tandem::order_and_save(int index[4], bool* state, int num_state1,
			int num_state2) {
    // First things first we need to order the 4 indexes i, j, k, l such
    // that i <= k and if i==k then j <= l by only doing changes of the
    // form i <-> j, k <-> l, (i, j) <-> (k, l).
    sort(index, index + 2);     // i <= j
    sort(index + 2, index + 4); // k <= l
    // demand i <= k. if not do (i, j) <-> (k, l).
    if (index[0] > index[2]) {
      swap(index[0], index[2]);
      swap(index[1], index[3]);
    }
    // i == j and j < l do (i, j) <-> (k, l).
    if ((index[0] == index[2]) && (index[1] > index[3])) {
      swap(index[0], index[2]);
      swap(index[1], index[3]);
    }
    // Now there is order, thats nice. But just to be sure.
    assert(index[0] <= index[2]);
    if (index[0] == index[2]) {
      assert(index[1] <= index[3]);
    }
    // Now we can sleep camly knowing nothing is weird with the
    // ordering.  Whats left is to determine the sign and then add
    // approriate wave_function elements to two_rdm.

    // The algoritm was designed for the erroneus formula.
    // p_ijkl = <phi|a^d_i*a^d_j*a_k*a_l|phi>
    // it turned out the correct formula is instead:
    // p_ijkl = <phi|a^d_i*a^d_j*a_l*a_k|phi>

    // Note the change k <-> l. This means this alogrithm is a (-1) sign
    // of... Instead of rewriting stuff and getting a lot of errors
    // I add a not infront of the sign statement.

    if (not sign(index, state)) {
      two_rdm[index[0]][index[1]][index[2]][index[3]] +=
          -wave_function[num_state1] * wave_function[num_state2];
    } else {
      two_rdm[index[0]][index[1]][index[2]][index[3]] +=
          wave_function[num_state1] * wave_function[num_state2];
    }
  }
