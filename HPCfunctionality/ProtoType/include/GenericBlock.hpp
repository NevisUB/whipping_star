#ifndef DIYGEN_GenericBlock_HPP
#define DIYGEN_GenericBlock_HPP

#include <math.h>

using namespace std;

// --------------- generic block definition ----------------
// B = Bounds definition
// S = State definition
// D = dataa definition

template <typename B, typename S, typename D>
struct GenericBlock
{
  typedef B bounds_type;
  typedef S state_type; 
  typedef D data_type;
  typedef GenericBlock<B,S,D> block_type;

  GenericBlock(bounds_type const& b):bounds(b),ah(NULL) { }

  // ---------------------------------------
  // standard functions for block processing
  static void* create() { return new block_type; }
  static void destroy(void* b) { delete static_cast<block_type*>(b); }

  static void save(const void* b, diy::BinaryBuffer& bb)
  {
    block_type const* bp = static_cast<block_type const*>(b);
    diy::save(bb, bp->bounds);
    diy::save(bb, bp->state);
    diy::save(bb, bp->data);
  }
  static void load(void* b, diy::BinaryBuffer& bb) 
  {
    block_type* bp = static_cast<block_type*>(b);
    diy::load(bb, bp->bounds);
    diy::load(bb, bp->state);
    diy::load(bb, bp->data);
  }

  // -------
  // this is the protocol for using the generic reduction function
  // this should be reduced to a minimum number of options

  // get access to the reduce buffer and the block state that is filled by processing
  data_type const& reduceData() const { return data; }
  data_type& reduceBuffer() { return buffer; }
  // a = a + b    (might want to permit other combining function to be passed in here)
  static void reduce(data_type& a, data_type const& b)
  {
    std::transform(a.cbegin(),a.cend(), b.cbegin(),
		   a.begin(), std::plus<typename data_type::value_type>());
  }
  // add "other" into my block data
  void reduce(data_type const& other) { reduce(data,other); }
  // add my block data into "other"
  void altreduce(data_type& other) const { reduce(other,data); }
  // add "other" into my buffer
  void reduce_buffer(data_type const& other) { reduce(buffer,other); }
  // add my buffer data into "other"
  void altreduce_buffer(data_type& other) const { altreduce(other,buffer); }

  // -----------
  // block data
  bounds_type bounds;
  state_type state;
  data_type data;
  data_type buffer;

private:
  GenericBlock() { }
};


#endif
