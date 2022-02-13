//
// Self-contained somewhat general-purpose hash table template class.
// The variable-length key must be a sequence of "atomic units" put in a std::vector container.
// The key "atomic unit" is specified by the 1st template parameter. 
// The value type associated with a key is specified by the 2nd template parameter.
//
// EXAMPLE:
//
//   // key = std::vector<int16_t>, value = int
//
//   hashtable<int16_t, int> HT(10);          // 2^10 = 1024 bucket empty hash table
//   hashkey<int16_t> k; k.push_back(12345);  // ... variable length key build up ...
//   int v = 1234;
//   HT.insert(k, v);                         // insert key-value pair (k, v)
//   bool hask = HT.exists(k);
//   HT.remove(k);
//
// Uses the MurmurHash3 hash function (which is public domain).
//

#ifndef __MUHASHTABLE_H__
#define __MUHASHTABLE_H__

#define __mumu3_hashsize(n) ((uint32_t)1 << (n))
#define __mumu3_hashmask(n) (__mumu3_hashsize(n) - 1)

static inline uint32_t __mumu3_rotl32(uint32_t x, int8_t r) {
  return (x << r) | (x >> (32 - r));
}

static inline uint32_t __mumu3_getblock32(const uint32_t* p, int i) {
  return p[i];
}

static inline uint32_t __mumu3_fmix32(uint32_t h) {
  h ^= h >> 16;
  h *= 0x85ebca6b;
  h ^= h >> 13;
  h *= 0xc2b2ae35;
  h ^= h >> 16;
  return h;
}

static inline uint32_t MurmurHash3_x86_32 (const void* key, int len, uint32_t seed) {
  const uint8_t* data = (const uint8_t*) key;
  const int nblocks = len / 4;
  uint32_t h1 = seed;
  const uint32_t c1 = 0xcc9e2d51;
  const uint32_t c2 = 0x1b873593;
  const uint32_t* blocks = (const uint32_t *)(data + nblocks * 4);
  for(int i = -nblocks; i; i++) {
    uint32_t k1 = __mumu3_getblock32(blocks, i);
    k1 *= c1;
    k1 = __mumu3_rotl32(k1, 15);
    k1 *= c2;
    h1 ^= k1;
    h1 = __mumu3_rotl32(h1, 13); 
    h1 = h1 * 5 + 0xe6546b64;
  }
  const uint8_t* tail = (const uint8_t*)(data + nblocks * 4);
  uint32_t k1 = 0;
  switch(len & 3)
  {
  case 3: k1 ^= tail[2] << 16;
  case 2: k1 ^= tail[1] << 8;
  case 1: k1 ^= tail[0];
          k1 *= c1; k1 = __mumu3_rotl32(k1, 15); k1 *= c2; h1 ^= k1;
  };
  h1 ^= len;
  h1 = __mumu3_fmix32(h1);
  return h1;
} 

template <typename A>
using hashkey = std::vector<A>;

template <typename A>
static inline bool __equal_keys(const hashkey<A>& k1, const hashkey<A>& k2) {
  if (k1.size() != k2.size()) return false;
  return (std::memcmp(k1.data(), k2.data(), sizeof(A) * k1.size()) == 0);
}

static inline int __upper_powtwo(int n, bool orequal = false) {
  int j = 0;
  if (orequal) while ((1 << j) < n) j++; // could stop so that 2^j == n
    else while ((1 << j) <= n) j++; // 2^j > n
  return j;
}

template<typename A, typename T>
class hashbucketitem
{
public:
  hashbucketitem(const hashkey<A>& k, const T& v) : key(k), value(v) { };
  hashkey<A> key;
  T value;
};

// A is the object that builds up the key in the std::vector container
// key is of type hashkey<A>, T is the value type
// the bucket item is a tuple {key, value} as defined above
template<typename A, typename T>
class hashtable
{
public:
  hashtable() { hashtable(10); };
  hashtable(int n, uint32_t seed = 0xdeadbeef) : bucket_(__mumu3_hashsize(n)) {
    items_ = 0;
    hashsize_ = n; // the bucket count will be 2^n
    seed_ = seed;
    bytes_per_atom_ = sizeof(A);
  };
  ~hashtable() {};

  int buckets() const { return static_cast<int>(bucket_.size()); };
  int items() const { return items_; };
  double load_factor() const { return static_cast<double>(items()) / buckets(); };
  bool empty() const { return (items_ == 0); };

  // return true if (key, value) did not exist and was inserted in its hash-bucket
  bool insert(const hashkey<A>& k, const T& v, bool front = true) {
    typename std::list<hashbucketitem<A, T> >* plb = &bucket_[hashbucket(k)];
    if (plb->empty()) {
      plb->push_front(hashbucketitem<A, T>(k, v));
      items_++;
      return true;
    }
    typename std::list<hashbucketitem<A, T> >::iterator it = plb->begin();
    while (it != plb->end()) {
      if (__equal_keys<A>(k, it->key)) return false;
      it++;
    }
    if (front) plb->push_front(hashbucketitem<A, T>(k, v));
      else plb->push_back(hashbucketitem<A, T>(k, v));
    items_++;
    return true;
  };

  // Locate key and remove the associated (k, v) if it is found
  bool remove(const hashkey<A>& k) {
    typename std::list<hashbucketitem<A, T> >* plb = &bucket_[hashbucket(k)];
    if (plb->empty()) return false;
    typename std::list<hashbucketitem<A, T> >::iterator it = plb->begin();
    while (it != plb->end()) {
      if (__equal_keys<A>(k, it->key)) {
        plb->erase(it);
        items_--;
        return true;
      }
      it++;
    }
    return false;
  };

  // Check whether the key k is stored in the table 
  bool exists(const hashkey<A>& k, int* pb = NULL) const {
    int bk = hashbucket(k);
    if (pb) *pb = bk;
    const typename std::list<hashbucketitem<A, T> >* plb = &bucket_[bk];
    if (plb->empty()) return false;
    typename std::list<hashbucketitem<A, T> >::const_iterator it = plb->begin();
    while (it != plb->end()) {
      if (__equal_keys<A>(k, it->key)) return true;
      it++;
    }
    return false;
  };

  // same as above but also extracts the value associated with the key (when returning true only)
  bool get(const hashkey<A>& k, T* v, int* pb = NULL) const {
    int bk = hashbucket(k);
    if (pb) *pb = bk;
    const typename std::list<hashbucketitem<A, T> >* plb = &bucket_[bk];
    if (plb->empty()) return false;
    typename std::list<hashbucketitem<A, T> >::const_iterator it = plb->begin();
    while (it != plb->end()) {
      if (__equal_keys<A>(k, it->key)) {
        *v = it->value;
        return true;
      }
      it++;
    }
    return false;
  };

  // pop() works like get(), but in addition pop() removes the matched key item
  bool pop(const hashkey<A>& k, T* v) {
    typename std::list<hashbucketitem<A, T> >* plb = &bucket_[hashbucket(k)];
    if (plb->empty()) return false;
    typename std::list<hashbucketitem<A, T> >::iterator it = plb->begin();
    while (it != plb->end()) {
      if (__equal_keys<A>(k, it->key)) {
        *v = it->value;
        plb->erase(it);
        items_--;
        return true;
      }
      it++;
    }
    return false;
  };

  void insert_front_nocheck(const hashkey<A>& k, const T& v) {
    insert_front_nocheck_bucket(hashbucket(k), k, v);
  };

  void insert_front_nocheck_bucket(int b, const hashkey<A>& k, const T& v) {
    bucket_[b].push_front(hashbucketitem<A, T>(k, v));
    items_++;
  };

  uint32_t hashcode(const hashkey<A>& k) const {
    return MurmurHash3_x86_32(reinterpret_cast<const void*>(k.data()), k.size() * bytes_per_atom_, seed_);
  };

  uint32_t hashbucket(const hashkey<A>& k) const {
    uint32_t h = hashcode(k);
    return (h & __mumu3_hashmask(hashsize_));
  }

  int bytes_per_atom() const { return bytes_per_atom_; };

  // Extract how many items are in each bucket
  int get_bucket_items(std::vector<int>* c, int* bmax) const {
    int csum = 0;
    if (bmax) *bmax = 0;
    if (c) c->clear();
    for (size_t i = 0; i < bucket_.size(); i++) {
      int bi = bucket_[i].size();
      if (c) c->push_back(bi);
      if (bmax) if (bi > *bmax) *bmax = bi;
      csum += bi;
    }
    return csum;
  };

  // Histogram such that v[i] = #buckets with exactly i items, i = 0...max < items()
  int get_bucket_histogram(std::vector<int>& h) const {
    std::vector<int> c;
    int bmax;
    int csum = get_bucket_items(&c, &bmax);
    h.clear();
    for (int i = 0; i <= bmax; i++) h.push_back(0);
    for (size_t i = 0; i < c.size(); i++) h[c[i]]++;
    return csum;
  };

  // Clear each bucket one by one; return the number of elements that were purged.
  int clear() {
    int csum = 0;
    for (size_t i = 0; i < bucket_.size(); i++) {
      if (!bucket_[i].empty()) {
        csum += bucket_[i].size();
        bucket_[i].clear();
      }
    }
    items_ -= csum; // should end up at zero exactly
    return csum;
  };

  // Basic resizing operations
  void upbuckets() {
    std::list<hashbucketitem<A, T> > L;
    exhale(L);
    hashsize_++;
    bucket_.resize(__mumu3_hashsize(hashsize_));
    inhale(L);
  };

  void downbuckets() {
    if (hashsize_ == 0) return; // only one bucket, minimally sized already
    std::list<hashbucketitem<A, T> > L;
    exhale(L);
    hashsize_--;
    bucket_.resize(__mumu3_hashsize(hashsize_));
    //bucket_.shrink_to_fit();
    inhale(L);
  };

  int rehash(uint32_t newseed) {
    if (items_ == 0) {
      seed_ = newseed;
      return items_;
    }
    std::list<hashbucketitem<A, T> > L;
    exhale(L);
    seed_ = newseed;
    inhale(L);
    return items_;
  };

  uint32_t seed() const { return seed_; };

  // transfer every item in the buckets onto a list (appending)
  void exhale(std::list<hashbucketitem<A, T> >& L) {
    for (size_t i = 0; i < bucket_.size(); i++) {
      typename std::list<hashbucketitem<A, T> >* pli = &bucket_[i];
      if (!pli->empty()) {
        typename std::list<hashbucketitem<A, T> >::iterator it = pli->begin();
        while (it != pli->end()) {
          L.push_back(hashbucketitem<A, T>(it->key, it->value));
          it = pli->erase(it);
          items_--; // final erase should result in 0 exactly
        }
      }
    }
  }

  // transfer every item in a list to the buckets
  void inhale(std::list<hashbucketitem<A, T> >& L) {
    typename std::list<hashbucketitem<A, T> >::iterator it = L.begin();
    while (it != L.end()) {
      bool was_inserted = insert(it->key, it->value);
      if (was_inserted) it = L.erase(it); else it++; // erase the accepted items only
    }
  }

  // Empty othertable items into this table with or without checking keys are unique
  void inhale(hashtable<A, T>& othertable, bool nocheck = false) {
    for (size_t i = 0; i < othertable.bucket_.size(); i++) {
      typename std::list<hashbucketitem<A, T> >* opli = &othertable.bucket_[i];
      if (!opli->empty()) {
        typename std::list<hashbucketitem<A, T> >::iterator oit = opli->begin();
        while (oit != opli->end()) {
          if (nocheck) {
            insert_front_nocheck(oit->key, oit->value); // no check that key is not in bucket already
            oit = opli->erase(oit);
            othertable.items_--;
            continue;
          }
          // only erase items from othertable that could be inserted into this table
          if (insert(oit->key, oit->value)) {
            oit = opli->erase(oit);
            othertable.items_--;
          } else {
            oit++;
          }
        }
      }
    }
  }

protected: 
  std::vector<std::list<hashbucketitem<A, T> > > bucket_; // forward_list would be enough really
  int hashsize_; // it will hold: bucket_.size() = 2^hashsize_
  int items_; // keep track of how many items are stored in the lists in all buckets
  uint32_t seed_; // which hash function ?
  int bytes_per_atom_;
};

#endif
