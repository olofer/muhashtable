//
// Basic test of Murmur3-based hashtable template class
//
// COMPILE: g++ -o test-muhashtable.exe -Wall -O2 test-muhashtable.cpp
//

#include <iostream>
#include <vector>
#include <list>
#include <cmath>   // for std::round(.)
#include <cstdlib> // for std::strtod(.)
#include <cstring> // for std::memcmp(.)

#include "muhashtable.h"

#include <random>
#include <chrono>

namespace RandUtils {
  std::mt19937 RandomGenerator;
  //std::default_random_engine RandomGenerator;
  void SetRandomGeneratorSeedExplicit(unsigned long ul) {
    RandomGenerator.seed(ul);
  }
  unsigned long SetRandomGeneratorSeedNow() {
    auto tp = std::chrono::high_resolution_clock::now();
    auto dn = tp.time_since_epoch();
    unsigned long ul = dn.count();
    RandomGenerator.seed(ul);
    return ul;
  }
  // Produce a single random unform integer in the given range [a, b]
  int UniformVariate(int a, int b) {
    std::uniform_real_distribution<double> U(0.0f, 1.0f);
    int j = std::round((b - a) * U(RandomGenerator)); // j is in [0, (b - a)]
    return (a + j);
  }
}

void generate_ideal_bucket_histogram(int N, int M, std::vector<int>& bh);

#define __HASHKEY_ATOM int32_t

int main(int argc, char** argv) 
{
  if (argc != 3) {
    std::cout << "usage: \"" << argv[0] << " n a\" (for 2^n buckets, load factor a)" << std::endl;
    return 0;
  }

  int n = static_cast<int>(std::round(std::strtod(argv[1], nullptr)));
  if (n <= 0) return 1;

  double alpha = std::strtod(argv[2], nullptr);
  if (alpha <= 0) return 1;

  unsigned long ul = RandUtils::SetRandomGeneratorSeedNow(); // try to pass this one thorugh as seed

  // hash table for "double" items
  // hashkeys are variable length vectors of __HASHKEY_ATOM's
  std::cout << "sizeof(__HASHKEY_ATOM) = " << sizeof(__HASHKEY_ATOM) << std::endl;

  std::vector<hashkey<__HASHKEY_ATOM>> keylog;

  hashtable<__HASHKEY_ATOM, double> HT(n, static_cast<uint32_t>(ul)); // ideally I want this template to take 2 parameters: <key unit type, value type>

  int m = static_cast<int>(std::round(alpha * HT.buckets()));

  std::cout << "Initially: HT.buckets() = " << HT.buckets() << ", HT.items() = " << HT.items() << std::endl;

  // go up to a load factor of 0.50
  for (int i = 0; i < m; i++) {
    hashkey<__HASHKEY_ATOM> ki;
    ki.push_back(0);
    ki.push_back(i); // this makes the complete key always unique [does not work great if the atom cannot handle large enough i:s]
    if (i & 0x01) ki.push_back(0); // every other key will have length 3 instead of 2 (# words)
    double vi = static_cast<double>(i);
    HT.insert(ki, vi);
    keylog.push_back(ki); // save keys so that I can remove everything
  }

  std::cout << "Post-insert: HT.buckets() = " << HT.buckets() << ", HT.items() = " << HT.items() << std::endl;

  std::vector<int> bh;
  HT.get_bucket_histogram(bh);
  for (size_t i = 0; i < bh.size(); i++) {
    std::cout << "# buckets with " << i << " items: " << bh[i] << std::endl;
  }

  std::cout << "Load factor = " << HT.load_factor() << std::endl;

  // key existence checks
  for (size_t i = 0; i < keylog.size(); i++) {
    if (!HT.exists(keylog[i])) {
      std::cout << "could not re-access key # " << i << "!" << std::endl;
      return 1;
    }
    hashkey<__HASHKEY_ATOM> ki;
    ki.push_back(i + 1);
    ki.push_back(i); // this key will never exist in HT
    if (HT.exists(ki)) {
      std::cout << "found key that was never inserted!" << std::endl;
      return 1;
    }
    // Now recover the value which should be equal to static_cast<double>(i) for keylog[i]
    double vi = -1.0;
    HT.get(keylog[i], &vi);
    if (vi != static_cast<double>(i)) {
      std::cout << "failed to recover the correct value for key # " << i << "!" << std::endl;
      return 1;
    }
  }
  std::cout << "passed key (non-)existence (and value recovery) checks" << std::endl;

  // Rehash with new seed
  unsigned long ul2 = RandUtils::SetRandomGeneratorSeedNow();
  HT.rehash(static_cast<uint32_t>(ul2));
  std::cout << "performed (constant-memory) rehash with new seed" << std::endl;

  // Check histogram again
  HT.get_bucket_histogram(bh);
  for (size_t i = 0; i < bh.size(); i++) {
    std::cout << "# buckets with " << i << " items: " << bh[i] << std::endl;
  }
  std::cout << "Load factor = " << HT.load_factor() << std::endl;

  // Resize experiments; upbucket(), show load factor, then downbucket(), show load factor
  HT.upbuckets();
  HT.get_bucket_histogram(bh);
  for (size_t i = 0; i < bh.size(); i++) {
    std::cout << "# buckets with " << i << " items: " << bh[i] << std::endl;
  }
  std::cout << "Load factor (post upbucket) = " << HT.load_factor() << std::endl;
  HT.downbuckets();
  HT.get_bucket_histogram(bh);
  for (size_t i = 0; i < bh.size(); i++) {
    std::cout << "# buckets with " << i << " items: " << bh[i] << std::endl;
  }
  std::cout << "Load factor (post downbucket) = " << HT.load_factor() << std::endl;

  // Next remove the (key,value) pairs by keys logged in keylog
  for (size_t i = 0; i < keylog.size(); i++) {
    if (!HT.remove(keylog[i])) {
      std::cout << "could not remove item with key # " << i << "!" << std::endl;
      return 1;
    }
  }

  std::cout << "Post-remove: HT.buckets() = " << HT.buckets() << ", HT.items() = " << HT.items() << std::endl;
  std::cout << "Load factor = " << HT.load_factor() << std::endl;

  HT.clear();

  std::cout << "*** comparison ideally uniform hashing sample ***" << std::endl;
  generate_ideal_bucket_histogram(HT.buckets(), m, bh);
  for (size_t i = 0; i < bh.size(); i++) {
    std::cout << "# buckets with " << i << " items: " << bh[i] << std::endl;
  }

  return 0;
}

// Subroutine for ideal uniform hashing comparison ...
void generate_ideal_bucket_histogram(int N, int M, std::vector<int>& bh) {
  std::vector<int> bin;
  for (int i = 0; i < N; i++) bin.push_back(0);
  int bmax = 0;
  for (int j = 0; j < M; j++) {
    int rj = RandUtils::UniformVariate(0, N - 1);
    bin[rj]++; // add one hit
    if (bin[rj] > bmax) bmax = bin[rj];
  }
  bh.clear();
  for (int i = 0; i <= bmax; i++) bh.push_back(0);
  for (int i = 0; i < N; i++) bh[bin[i]]++;
  return;
}
