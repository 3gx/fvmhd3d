#ifndef _MEMORY_POOL_H_
#define _MEMORY_POOL_H_

#include <vector>

template<typename T>
inline bool safe_resize(std::vector<T> &vec, size_t size) {
	T* ptr_old = &vec[0];
//	erase_vec(vec);
	vec.resize(size);
	return (ptr_old == &vec[0]);
}

template<class T>
struct memory_pool {


  std::vector<T> vec;
  size_t capacity;
  size_t used;


  memory_pool() {
    vec.clear();
    capacity = 0;
    used = 0;
  };
  ~memory_pool() {};
  
  /////////
  
  
  bool resize(const size_t new_capacity) {
#if 0
    if (new_capacity > capacity) {
      capacity = new_capacity;
      return safe_resize(vec, capacity);
    }
#else
		capacity = new_capacity;
		return safe_resize(vec, capacity);
#endif
    return true;
  }

	bool free_unused() {
		return true;
		T* ptr = &vec[0];
		vec.resize(used);
//		vec.reserve(used);
		return ptr == &vec[0];
	}

  T *allocate(size_t ss) {
    if (ss > capacity) {
      vec.resize(ss);
//      vec.reserve(ss);
      capacity = ss;
    }
    used = 0;
    return &vec[0];
  }

  T *get(const int count) {
    T *ret = &vec[used];
    used += count;
    if (!(used <= capacity)) {
      fprintf(stderr,  " used= %d > capacity= %d \n", (int)used, (int)capacity);
    }
    assert(used <= capacity);
    
    return ret;
  };

  T& operator[](int i){return vec[i];}
  void reset() {
    used = 0;
  }

};


#endif // _MEMORY_POOL_H_
