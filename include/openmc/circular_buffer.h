#ifndef OPENMC_CIRCULAR_BUFFER_H
#define OPENMC_CIRCULAR_BUFFER_H

#include <algorithm>
#include <stdexcept>
#include <vector>

namespace openmc {

//! A fixed-size circular buffer that provides vector-like interface
//! with automatic handling of pushing new elements.
template<typename T>
class CircularBuffer {
public:
  //! Default constructor - creates an empty buffer
  CircularBuffer() : data_(), capacity_(0), size_(0), start_(0) {}

  //! Constructor with specified capacity
  explicit CircularBuffer(size_t capacity)
    : data_(capacity), capacity_(capacity), size_(0), start_(0)
  {}

  //! Copy constructor
  CircularBuffer(const CircularBuffer& other) = default;

  //! Constructor from vector - initializes with vector content
  explicit CircularBuffer(const std::vector<T>& vec)
    : data_(vec), capacity_(vec.size()), size_(vec.size()), start_(0)
  {}

  //! Assignment operator from vector
  CircularBuffer& operator=(const std::vector<T>& vec)
  {
    // Clear existing data
    clear();

    // If vector is larger than our capacity, resize
    if (vec.size() > capacity_) {
      data_.resize(vec.size());
      capacity_ = vec.size();
    }

    // Copy elements from vector
    for (const auto& value : vec) {
      push_back(value);
    }

    return *this;
  }

  //! Access element with circular indexing
  T& operator[](size_t index)
  {
    if (index >= size_)
      throw std::out_of_range("CircularBuffer index out of range");
    return data_[(start_ + index) % capacity_];
  }

  //! Const access element with circular indexing
  const T& operator[](size_t index) const
  {
    if (index >= size_)
      throw std::out_of_range("CircularBuffer index out of range");
    return data_[(start_ + index) % capacity_];
  }

  //! Return the current size
  size_t size() const { return size_; }

  //! Return the maximum capacity
  size_t capacity() const { return capacity_; }

  //! Check if buffer is empty
  bool empty() const { return size_ == 0; }

  //! Resize the buffer (similar to vector's resize)
  void resize(size_t n)
  {
    if (n <= capacity_) {
      size_ = n;
      start_ = 0; // Simplify by resetting start position
    } else {
      // Need to increase capacity
      data_.resize(n);
      capacity_ = n;
      size_ = n;
      start_ = 0;
    }
  }

  //! Add an element to the end, handle circular behavior
  void push_back(const T& value)
  {
    if (size_ < capacity_) {
      // Still have room in buffer
      data_[(start_ + size_) % capacity_] = value;
      ++size_;
    } else {
      // Buffer is full, remove first element and add new one
      data_[start_] = value;
      start_ = (start_ + 1) % capacity_;
    }
  }

  //! Remove the oldest element
  void pop_front()
  {
    if (size_ > 0) {
      start_ = (start_ + 1) % capacity_;
      --size_;
    }
  }

  //! Clear the buffer
  void clear()
  {
    size_ = 0;
    start_ = 0;
  }

  //! Convert to vector (useful for compatibility with existing code)
  std::vector<T> to_vector() const
  {
    std::vector<T> result(size_);
    for (size_t i = 0; i < size_; ++i) {
      result[i] = data_[(start_ + i) % capacity_];
    }
    return result;
  }

private:
  std::vector<T> data_; // Underlying storage
  size_t capacity_;     // Maximum capacity
  size_t size_;         // Current number of elements
  size_t start_;        // Index of the first element
};

} // namespace openmc

#endif // OPENMC_CIRCULAR_BUFFER_H
