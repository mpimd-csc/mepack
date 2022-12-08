#pragma once

#include "matrix.h"
#include <limits>

// A wrapper for matlabs Allocator to allow usage of std containers.
template<typename T>
class MexAllocator {
public : 
	typedef T value_type;
	typedef value_type* pointer;
	typedef const value_type* const_pointer;
	typedef value_type& reference;
	typedef const value_type& const_reference;
	typedef std::size_t size_type;
	typedef std::ptrdiff_t difference_type;

public : 
	//    convert an allocator<T> to allocator<U>
	template<typename U>
	struct rebind {
		typedef MexAllocator<U> other;
	};

public : 
	explicit MexAllocator() {}
	~MexAllocator() {}
	MexAllocator(MexAllocator const&) {}
	template<typename U>
	explicit MexAllocator(MexAllocator<U> const&) {}

	pointer allocate(size_type cnt, typename std::allocator<void>::const_pointer = 0) 
	{ 
	  return reinterpret_cast<pointer>(mxMalloc(cnt * sizeof (T))); 
	}
	
	void deallocate(pointer p, size_type) 
	{ 
		mxFree(p); 
	}

	inline size_type max_size() const 
	{ 
		return std::numeric_limits<size_type>::max() / sizeof(T);
 	}

	void construct(pointer p, const T& t) { new(p) T(t); }
	void destroy(pointer p) { p->~T(); }

	inline bool operator==(MexAllocator const&) { return true; }
	inline bool operator!=(MexAllocator const& a) { return !operator==(a); }
}; 
