#pragma once
#include <array>
#include <cassert>
#include <complex>
#include <initializer_list>
#include <memory>
#include <vector>

#ifndef __OPTIMIZE__
  #define kBoundsCheck
#endif

static const int kMaxArrayDimension = 11;

template <typename T>
struct TArray {
  int nel = 0;
  int ndim = 0;
  std::array<int, kMaxArrayDimension> dims;
  std::array<int, kMaxArrayDimension> strides;
  bool ownsData = false;
  T *data = nullptr;

  TArray() = default;
  TArray(const std::vector<int> &ns) { init(ns); alloc(); }
  TArray(T *ptr, const std::vector<int> &ns) {
    init(ns);
    data = ptr;
  }
  TArray(std::nullptr_t, const std::vector<int> &ns) { init(ns); }
  template <typename... D>
  TArray(std::nullptr_t, D... ns) {
    init(ns...);
  }
  template <typename... D>
  TArray(T *ptr, D... ns) {
    init(ns...);
    data = ptr;
  }
  template <typename... D>
  TArray(D... ns) {
    init(ns...);
    alloc();
  }
  template <typename ...D>
  void init(T *ptr, D... ns) {
    init(ns...);
    data = ptr;
  }
  template <typename ...D>
  void init(D... ns) {
    dims = {{static_cast<int>(ns)...}};
    ndim = sizeof...(ns);
    init();
  }
  void init(const std::vector<int> &ns) {
    ndim = ns.size();
    for (int i=0; i<kMaxArrayDimension; ++i)
      dims[i] = 0;
    std::copy(ns.begin(), ns.end(), dims.begin());
    init();
  }
  void init() {
    nel = 1;
    strides[0] = 1;
    for (int i=0; i<kMaxArrayDimension; ++i) {
      int d = dims[i];
      if (d == 0) break;
      nel *= d;
      strides[i+1] = nel;
    }
  }

  template <typename ...D>
  void alloc(D... ns) {
    init(ns...);
    alloc();
  }
  void alloc(const std::vector<int> &ns) {
    init(ns);
    alloc();
  }
  void alloc() {
    #ifdef kBoundsCheck
      assert(!ownsData);
    #endif
    ownsData = true;
    data = new T[nel];
  }

  void realloc(const std::vector<int> &newDims) {
    if (!ownsData) {
      alloc(newDims);
    } else {
      bool same = true;
      if (size_t(ndim) == newDims.size()) {
        for (int d=0; d<ndim; ++d) {
          if (dims[d] != newDims[d]) {
            same = false;
            break;
          }
        }
      } else {
        same = false;
      }
      if (!same) {
        clear();
        alloc(newDims);
      }
    }
  }

  void setData(T *ptr) {
    #ifdef kBoundsCheck
      assert(!ownsData);
    #endif
    data = ptr;
  }

  std::vector<int> size() const {
    std::vector<int> result(ndim);
    std::copy(dims.begin(), dims.begin()+ndim, result.begin());
    return result;
  }
  int size(int d) const {
    #ifdef kBoundsCheck
      assert(d < ndim);
    #endif
    return dims[d];
  }

  operator T*() const { return data; }
  T* begin() const { return data; }
  T* end() const { return data + nel; }

  T &operator[](int i) const {
    #ifdef kBoundsCheck
      assert(i < nel);
    #endif
    return data[i*strides[0]];
  }

  T &index(int i, int d) const {
    #ifdef kBoundsCheck
      assert(d == ndim);
    #endif
    return data[i];
  }
  template <typename... Idx>
  T &index(int i, int d, int first, Idx... rest) const {
    #ifdef kBoundsCheck
      assert(d < ndim);
      assert(first < dims[d] && first >= 0);
    #endif
    return index(i + first*strides[d], d+1, rest...);
  }
  template <typename... Idx>
  T &operator()(int first, Idx... rest) const {
    #ifdef kBoundsCheck
      assert(first < dims[0] && first >= 0);
    #endif
    return index(first*strides[0], 1, rest...);
  }

  void operator=(T x) {
    for (int i=0; i<nel; ++i) {
      data[i*strides[0]] = x;
    }
  }

  void clear() {
    nel = 0;
    ndim = 0;
    dims.fill(0);
    strides.fill(0);
    if (ownsData && data) {
      delete [] data;
    }
    data = nullptr;
    ownsData = false;
  }



  ~TArray() {
    if (ownsData && data) delete [] data;
  }

  //////////////////////////////////////////////////////////////////////////////
  // Array creation using initializer lists (a bit complicated...)            //
  //////////////////////////////////////////////////////////////////////////////
  struct ArrayInitList {
    std::shared_ptr<void> subLists; // empty => d is valid
    T d;
    ArrayInitList(T v) : d(v) { }
    ArrayInitList(std::initializer_list<ArrayInitList> l)
     : subLists(std::make_shared<std::vector<ArrayInitList>>(l)), d() {
    }
  };

  TArray& operator=(std::initializer_list<ArrayInitList> l) {
    clear();
    processInitList(ArrayInitList(l));
    return *this;
  }

  TArray(std::initializer_list<ArrayInitList> l) {
    processInitList(ArrayInitList(l));
  }

  void processInitList(ArrayInitList l, std::vector<int> idx = {}) {
    if (l.subLists) {
      std::shared_ptr<std::vector<ArrayInitList>> subLists =
        std::static_pointer_cast<std::vector<ArrayInitList>>(l.subLists);
      if (idx.size() >= ndim) {
        ++ndim;
        strides[ndim-1] = 1;
        dims[ndim-1] = subLists->size();
        for (int j=0; j<ndim-1; ++j)
          strides[ndim-1] *= dims[j];
      } else {
        assert(dims[idx.size()] == subLists->size());
      }
      idx.push_back(0);
      for (const auto &sl : *subLists) {
        processInitList(sl, idx);
        ++idx.back();
      }
    } else {
      if (!data) {
        ownsData = true;
        nel = 1;
        for (int j=0; j<ndim; ++j) nel *= dims[j];
        data = new T[nel];
      }
      int linearIdx = 0;
      for (int j=0; j<ndim; ++j) linearIdx += idx[j]*strides[j];
      data[linearIdx] = l.d;
    }
  }
  //////////////////////////////////////////////////////////////////////////////
  // End of array creation using initializer lists                            //
  //////////////////////////////////////////////////////////////////////////////

  // Arrays are movable
  TArray(TArray &&x) noexcept
   : nel(x.nel), ndim(x.ndim), dims(std::move(x.dims)),
     strides(std::move(x.strides)), ownsData(x.ownsData), data(x.data) {
    x.ownsData = false;
    x.data = nullptr;
    x.nel = 0;
    x.ndim = 0;
  }
  TArray& operator=(TArray &&x) noexcept {
    nel = x.nel;
    ndim = x.ndim;
    dims = std::move(x.dims);
    strides = std::move(x.strides);
    ownsData = x.ownsData;
    data = x.data;
    x.ownsData = false;
    x.data = nullptr;
    x.nel = 0;
    x.ndim = 0;
    return *this;
  }
  // Prohibit copy and assignment operators
  TArray(TArray &) = delete;
  TArray& operator=(TArray &) = delete;
};

using DArray = TArray<double>;
using IArray = TArray<int>;
using ZArray = TArray<std::complex<double>>;

constexpr int size_product(int first) {
  return first;
}

template<typename... D>
constexpr int size_product(int first, D... rest) {
  return first*size_product(rest...);
}

template <typename T, int... D>
struct TArrayFixed : TArray<T> {
  T fixedData[size_product(D...)];
  TArrayFixed() : TArray<T>(nullptr, D...) {
    this->data = fixedData;
  }
  using TArray<T>::operator=;
};

template <int... D>
using DArrayFixed = TArrayFixed<double, D...>;
template <int... D>
using IArrayFixed = TArrayFixed<int, D...>;
template <int... D>
using ZArrayFixed = TArrayFixed<std::complex<double>, D...>;
