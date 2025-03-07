#if !defined(sparsepp_h_guard_)
#define sparsepp_h_guard_




// some macros for portability
// ---------------------------
// includes
// --------
#include <cassert>
#include <cstring>
#include <string>
#include <limits>                           // for numeric_limits
#include <algorithm>                        // For swap(), eg
#include <iterator>                         // for iterator tags
#include <functional>                       // for equal_to<>, select1st<>, std::unary_function, etc
#include <memory>                           // for alloc, uninitialized_copy, uninitialized_fill
#include <cstdlib>                          // for malloc/realloc/free
#include <cstddef>                          // for ptrdiff_t
#include <new>                              // for placement new
#include <stdexcept>                        // For length_error
#include <utility>                          // for pair<>
#include <cstdio>
#include <iosfwd>
#include <ios>

#include "spp_stdint.h"  // includes spp_config.h
#include "spp_traits.h"
#include "spp_utils.h"

#ifdef SPP_INCLUDE_SPP_ALLOC
    #include "spp_dlalloc.h"
#endif

#if !defined(SPP_NO_CXX11_HDR_INITIALIZER_LIST)
    #include <initializer_list>
#endif

#if (SPP_GROUP_SIZE == 32)
    #define SPP_SHIFT_ 5
    #define SPP_MASK_  0x1F
    typedef uint32_t group_bm_type;
#elif (SPP_GROUP_SIZE == 64)
    #define SPP_SHIFT_ 6
    #define SPP_MASK_  0x3F
    typedef uint64_t group_bm_type;
#else
    #error "SPP_GROUP_SIZE must be either 32 or 64"
#endif

namespace spp_ {

//  ----------------------------------------------------------------------
//                  U T I L    F U N C T I O N S
//  ----------------------------------------------------------------------
template <class E>
inline void throw_exception(const E& exception)
{
#if !defined(SPP_NO_EXCEPTIONS)
    throw exception;
#else
    assert(0);
    abort();
#endif
}

//  ----------------------------------------------------------------------
//              M U T A B L E     P A I R      H A C K
// turn std::pair<const K, V> into mutable std::pair<K, V>
//  ----------------------------------------------------------------------
template <class T>
struct cvt
{
    typedef T type;
};

template <class K, class V>
struct cvt<std::pair<const K, V> >
{
    typedef std::pair<K, V> type;
};

template <class K, class V>
struct cvt<const std::pair<const K, V> >
{
    typedef const std::pair<K, V> type;
};

//  ----------------------------------------------------------------------
//              M O V E   I T E R A T O R
//  ----------------------------------------------------------------------
#ifdef SPP_NO_CXX11_RVALUE_REFERENCES
    #define MK_MOVE_IT(p) (p)
#else
    #define MK_MOVE_IT(p) std::make_move_iterator(p)
#endif


//  ----------------------------------------------------------------------
//             I N T E R N A L    S T U F F
//  ----------------------------------------------------------------------
#ifdef SPP_NO_CXX11_STATIC_ASSERT
    template <bool> struct SppCompileAssert { };
    #define SPP_COMPILE_ASSERT(expr, msg) \
      SPP_ATTRIBUTE_UNUSED typedef SppCompileAssert<(bool(expr))> spp_bogus_[bool(expr) ? 1 : -1]
#else
    #define SPP_COMPILE_ASSERT static_assert
#endif

namespace sparsehash_internal
{

// Adaptor methods for reading/writing data from an INPUT or OUPTUT
// variable passed to serialize() or unserialize().  For now we
// have implemented INPUT/OUTPUT for FILE*, istream*/ostream* (note
// they are pointers, unlike typical use), or else a pointer to
// something that supports a Read()/Write() method.
//
// For technical reasons, we implement read_data/write_data in two
// stages.  The actual work is done in *_data_internal, which takes
// the stream argument twice: once as a template type, and once with
// normal type information.  (We only use the second version.)  We do
// this because of how C++ picks what function overload to use.  If we
// implemented this the naive way:
//    bool read_data(istream* is, const void* data, size_t length);
//    template<typename T> read_data(T* fp,  const void* data, size_t length);
// C++ would prefer the second version for every stream type except
// istream.  However, we want C++ to prefer the first version for
// streams that are *subclasses* of istream, such as istringstream.
// This is not possible given the way template types are resolved.  So
// we split the stream argument in two, one of which is templated and
// one of which is not.  The specialized functions (like the istream
// version above) ignore the template arg and use the second, 'type'
// arg, getting subclass matching as normal.  The 'catch-all'
// functions (the second version above) use the template arg to deduce
// the type, and use a second, void* arg to achieve the desired
// 'catch-all' semantics.

    // ----- low-level I/O for FILE* ----

    template<typename Ignored>
    inline bool read_data_internal(Ignored* /*unused*/, FILE* fp,
                                   void* data, size_t length)
    {
        return fread(data, length, 1, fp) == 1;
    }

    template<typename Ignored>
    inline bool write_data_internal(Ignored* /*unused*/, FILE* fp,
                                    const void* data, size_t length)
    {
        return fwrite(data, length, 1, fp) == 1;
    }

    // ----- low-level I/O for iostream ----

    // We want the caller to be responsible for #including <iostream>, not
    // us, because iostream is a big header!  According to the standard,
    // it's only legal to delay the instantiation the way we want to if
    // the istream/ostream is a template type.  So we jump through hoops.
    template<typename ISTREAM>
    inline bool read_data_internal_for_istream(ISTREAM* fp,
                                               void* data, size_t length)
    {
        return fp->read(reinterpret_cast<char*>(data),
                        static_cast<std::streamsize>(length)).good();
    }
    template<typename Ignored>
    inline bool read_data_internal(Ignored* /*unused*/, std::istream* fp,
                                   void* data, size_t length)
    {
        return read_data_internal_for_istream(fp, data, length);
    }

    template<typename OSTREAM>
    inline bool write_data_internal_for_ostream(OSTREAM* fp,
                                                const void* data, size_t length)
    {
        return fp->write(reinterpret_cast<const char*>(data),
                         static_cast<std::streamsize>(length)).good();
    }
    template<typename Ignored>
    inline bool write_data_internal(Ignored* /*unused*/, std::ostream* fp,
                                    const void* data, size_t length)
    {
        return write_data_internal_for_ostream(fp, data, length);
    }

    // ----- low-level I/O for custom streams ----

    // The INPUT type needs to support a Read() method that takes a
    // buffer and a length and returns the number of bytes read.
    template <typename INPUT>
    inline bool read_data_internal(INPUT* fp, void* /*unused*/,
                                   void* data, size_t length)
    {
        return static_cast<size_t>(fp->Read(data, length)) == length;
    }

    // The OUTPUT type needs to support a Write() operation that takes
    // a buffer and a length and returns the number of bytes written.
    template <typename OUTPUT>
    inline bool write_data_internal(OUTPUT* fp, void* /*unused*/,
                                    const void* data, size_t length)
    {
        return static_cast<size_t>(fp->Write(data, length)) == length;
    }

    // ----- low-level I/O: the public API ----

    template <typename INPUT>
    inline bool read_data(INPUT* fp, void* data, size_t length)
    {
        return read_data_internal(fp, fp, data, length);
    }

    template <typename OUTPUT>
    inline bool write_data(OUTPUT* fp, const void* data, size_t length)
    {
        return write_data_internal(fp, fp, data, length);
    }

    // Uses read_data() and write_data() to read/write an integer.
    // length is the number of bytes to read/write (which may differ
    // from sizeof(IntType), allowing us to save on a 32-bit system
    // and load on a 64-bit system).  Excess bytes are taken to be 0.
    // INPUT and OUTPUT must match legal inputs to read/write_data (above).
    // --------------------------------------------------------------------
    template <typename INPUT, typename IntType>
    bool read_bigendian_number(INPUT* fp, IntType* value, size_t length)
    {
        *value = 0;
        unsigned char byte;
        // We require IntType to be unsigned or else the shifting gets all screwy.
        SPP_COMPILE_ASSERT(static_cast<IntType>(-1) > static_cast<IntType>(0), "serializing_int_requires_an_unsigned_type");
        for (size_t i = 0; i < length; ++i)
        {
            if (!read_data(fp, &byte, sizeof(byte)))
                return false;
            *value |= static_cast<IntType>(byte) << ((length - 1 - i) * 8);
        }
        return true;
    }

    template <typename OUTPUT, typename IntType>
    bool write_bigendian_number(OUTPUT* fp, IntType value, size_t length)
    {
        unsigned char byte;
        // We require IntType to be unsigned or else the shifting gets all screwy.
        SPP_COMPILE_ASSERT(static_cast<IntType>(-1) > static_cast<IntType>(0), "serializing_int_requires_an_unsigned_type");
        for (size_t i = 0; i < length; ++i)
        {
            byte = (sizeof(value) <= length-1 - i)
                ? static_cast<unsigned char>(0) : static_cast<unsigned char>((value >> ((length-1 - i) * 8)) & 255);
            if (!write_data(fp, &byte, sizeof(byte))) return false;
        }
        return true;
    }

    // If your keys and values are simple enough, you can pass this
    // serializer to serialize()/unserialize().  "Simple enough" means
    // value_type is a POD type that contains no pointers.  Note,
    // however, we don't try to normalize endianness.
    // This is the type used for NopointerSerializer.
    // ---------------------------------------------------------------
    template <typename value_type> struct pod_serializer
    {
        template <typename INPUT>
        bool operator()(INPUT* fp, value_type* value) const
        {
            return read_data(fp, value, sizeof(*value));
        }

        template <typename OUTPUT>
        bool operator()(OUTPUT* fp, const value_type& value) const
        {
            return write_data(fp, &value, sizeof(value));
        }
    };


    // Settings contains parameters for growing and shrinking the table.
    // It also packages zero-size functor (ie. hasher).
    //
    // It does some munging of the hash value for the cases where
    // the original hash function is not be very good.
    // ---------------------------------------------------------------
    template<typename Key, typename HashFunc, typename SizeType, int HT_MIN_BUCKETS>
    class sh_hashtable_settings : public HashFunc
    {
    private:
#ifndef SPP_MIX_HASH
        template <class T, int sz> struct Mixer
        {
            inline T operator()(T h) const { return h; }
        };
#else
        template <class T, int sz> struct Mixer
        {
            inline T operator()(T h) const;
        };

         template <class T> struct Mixer<T, 4>
        {
            inline T operator()(T h) const
            {
                // from Thomas Wang - https://gist.github.com/badboy/6267743
                // ---------------------------------------------------------
                h = (h ^ 61) ^ (h >> 16);
                h = h + (h << 3);
                h = h ^ (h >> 4);
                h = h * 0x27d4eb2d;
                h = h ^ (h >> 15);
                return h;
            }
        };

        template <class T> struct Mixer<T, 8>
        {
            inline T operator()(T h) const
            {
                // from Thomas Wang - https://gist.github.com/badboy/6267743
                // ---------------------------------------------------------
                h = (~h) + (h << 21);              // h = (h << 21) - h - 1;
                h = h ^ (h >> 24);
                h = (h + (h << 3)) + (h << 8);     // h * 265
                h = h ^ (h >> 14);
                h = (h + (h << 2)) + (h << 4);     // h * 21
                h = h ^ (h >> 28);
                h = h + (h << 31);
                return h;
            }
        };
#endif

    public:
        typedef Key key_type;
        typedef HashFunc hasher;
        typedef SizeType size_type;

    public:
        sh_hashtable_settings(const hasher& hf,
                              const float ht_occupancy_flt,
                              const float ht_empty_flt)
            : hasher(hf),
              enlarge_threshold_(0),
              shrink_threshold_(0),
              consider_shrink_(false),
              num_ht_copies_(0)
        {
            set_enlarge_factor(ht_occupancy_flt);
            set_shrink_factor(ht_empty_flt);
        }

        size_t hash(const key_type& v) const
        {
            size_t h = hasher::operator()(v);
            Mixer<size_t, sizeof(size_t)> mixer;

            return mixer(h);
        }

        float enlarge_factor() const            { return enlarge_factor_; }
        void set_enlarge_factor(float f)        { enlarge_factor_ = f;    }
        float shrink_factor() const             { return shrink_factor_;  }
        void set_shrink_factor(float f)         { shrink_factor_ = f;     }

        size_type enlarge_threshold() const     { return enlarge_threshold_; }
        void set_enlarge_threshold(size_type t) { enlarge_threshold_ = t; }
        size_type shrink_threshold() const      { return shrink_threshold_; }
        void set_shrink_threshold(size_type t)  { shrink_threshold_ = t; }

        size_type enlarge_size(size_type x) const { return static_cast<size_type>(x * enlarge_factor_); }
        size_type shrink_size(size_type x) const { return static_cast<size_type>(x * shrink_factor_); }

        bool consider_shrink() const            { return consider_shrink_; }
        void set_consider_shrink(bool t)        { consider_shrink_ = t; }

        unsigned int num_ht_copies() const      { return num_ht_copies_; }
        void inc_num_ht_copies()                { ++num_ht_copies_; }

        // Reset the enlarge and shrink thresholds
        void reset_thresholds(size_type num_buckets)
        {
            set_enlarge_threshold(enlarge_size(num_buckets));
            set_shrink_threshold(shrink_size(num_buckets));
            // whatever caused us to reset already considered
            set_consider_shrink(false);
        }

        // Caller is resposible for calling reset_threshold right after
        // set_resizing_parameters.
        // ------------------------------------------------------------
        void set_resizing_parameters(float shrink, float grow)
        {
            assert(shrink >= 0);
            assert(grow <= 1);
            if (shrink > grow/2.0f)
                shrink = grow / 2.0f;     // otherwise we thrash hashtable size
            set_shrink_factor(shrink);
            set_enlarge_factor(grow);
        }

        // This is the smallest size a hashtable can be without being too crowded
        // If you like, you can give a min #buckets as well as a min #elts
        // ----------------------------------------------------------------------
        size_type min_buckets(size_type num_elts, size_type min_buckets_wanted)
        {
            float enlarge = enlarge_factor();
            size_type sz = HT_MIN_BUCKETS;             // min buckets allowed
            while (sz < min_buckets_wanted ||
                   num_elts >= static_cast<size_type>(sz * enlarge))
            {
                // This just prevents overflowing size_type, since sz can exceed
                // max_size() here.
                // -------------------------------------------------------------
                if (static_cast<size_type>(sz * 2) < sz)
                    throw_exception(std::length_error("resize overflow"));  // protect against overflow
                sz *= 2;
            }
            return sz;
        }

    private:
        size_type enlarge_threshold_;  // table.size() * enlarge_factor
        size_type shrink_threshold_;   // table.size() * shrink_factor
        float enlarge_factor_;         // how full before resize
        float shrink_factor_;          // how empty before resize
        bool consider_shrink_;         // if we should try to shrink before next insert

        unsigned int num_ht_copies_;   // num_ht_copies is a counter incremented every Copy/Move
    };

}  // namespace sparsehash_internal

#undef SPP_COMPILE_ASSERT

//  ----------------------------------------------------------------------
//                    S P A R S E T A B L E
//  ----------------------------------------------------------------------
//
// A sparsetable is a random container that implements a sparse array,
// that is, an array that uses very little memory to store unassigned
// indices (in this case, between 1-2 bits per unassigned index).  For
// instance, if you allocate an array of size 5 and assign a[2] = <big
// struct>, then a[2] will take up a lot of memory but a[0], a[1],
// a[3], and a[4] will not.  Array elements that have a value are
// called "assigned".  Array elements that have no value yet, or have
// had their value cleared using erase() or clear(), are called
// "unassigned".
//
// Unassigned values seem to have the default value of T (see below).
// Nevertheless, there is a difference between an unassigned index and
// one explicitly assigned the value of T().  The latter is considered
// assigned.
//
// Access to an array element is constant time, as is insertion and
// deletion.  Insertion and deletion may be fairly slow, however:
// because of this container's memory economy, each insert and delete
// causes a memory reallocation.
//
// NOTE: You should not test(), get(), or set() any index that is
// greater than sparsetable.size().  If you need to do that, call
// resize() first.
//
// --- Template parameters
// PARAMETER   DESCRIPTION                           DEFAULT
// T           The value of the array: the type of   --
//             object that is stored in the array.
//
// Alloc:      Allocator to use to allocate memory.
//
// --- Model of
// Random Access Container
//
// --- Type requirements
// T must be Copy Constructible. It need not be Assignable.
//
// --- Public base classes
// None.
//
// --- Members
//
// [*] All iterators are const in a sparsetable (though nonempty_iterators
//     may not be).  Use get() and set() to assign values, not iterators.
//
// [+] iterators are random-access iterators.  nonempty_iterators are
//     bidirectional iterators.

// [*] If you shrink a sparsetable using resize(), assigned elements
// past the end of the table are removed using erase().  If you grow
// a sparsetable, new unassigned indices are created.
//
// [+] Note that operator[] returns a const reference.  You must use
// set() to change the value of a table element.
//
// [!] Unassignment also calls the destructor.
//
// Iterators are invalidated whenever an item is inserted or
// deleted (ie set() or erase() is used) or when the size of
// the table changes (ie resize() or clear() is used).



// ---------------------------------------------------------------------------
// Our iterator as simple as iterators can be: basically it's just
// the index into our table.  Dereference, the only complicated
// thing, we punt to the table class.  This just goes to show how
// much machinery STL requires to do even the most trivial tasks.
//
// A NOTE ON ASSIGNING:
// A sparse table does not actually allocate memory for entries
// that are not filled.  Because of this, it becomes complicated
// to have a non-const iterator: we don't know, if the iterator points
// to a not-filled bucket, whether you plan to fill it with something
// or whether you plan to read its value (in which case you'll get
// the default bucket value).  Therefore, while we can define const
// operations in a pretty 'normal' way, for non-const operations, we
// define something that returns a helper object with operator= and
// operator& that allocate a bucket lazily.  We use this for table[]
// and also for regular table iterators.

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
// Our iterator as simple as iterators can be: basically it's just
// the index into our table.  Dereference, the only complicated
// thing, we punt to the table class.  This just goes to show how
// much machinery STL requires to do even the most trivial tasks.
//
// By templatizing over tabletype, we have one iterator type which
// we can use for both sparsetables and sparsebins.  In fact it
// works on any class that allows size() and operator[] (eg vector),
// as long as it does the standard STL typedefs too (eg value_type).

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
template <class tabletype>
class table_iterator
{
public:
    typedef table_iterator iterator;

    typedef std::random_access_iterator_tag      iterator_category;
    typedef typename tabletype::value_type       value_type;
    typedef typename tabletype::difference_type  difference_type;
    typedef typename tabletype::size_type        size_type;

    explicit table_iterator(tabletype *tbl = 0, size_type p = 0) :
        table(tbl), pos(p)
    { }

    // Helper function to assert things are ok; eg pos is still in range
    void check() const
    {
        assert(table);
        assert(pos <= table->size());
    }

    // Arithmetic: we just do arithmetic on pos.  We don't even need to
    // do bounds checking, since STL doesn't consider that its job.  :-)
    iterator& operator+=(size_type t) { pos += t; check(); return *this; }
    iterator& operator-=(size_type t) { pos -= t; check(); return *this; }
    iterator& operator++()            { ++pos; check(); return *this; }
    iterator& operator--()            { --pos; check(); return *this; }
    iterator operator++(int)
    {
        iterator tmp(*this);     // for x++
        ++pos; check(); return tmp;
    }

    iterator operator--(int)
    {
        iterator tmp(*this);     // for x--
        --pos; check(); return tmp;
    }

    iterator operator+(difference_type i) const
    {
        iterator tmp(*this);
        tmp += i; return tmp;
    }

    iterator operator-(difference_type i) const
    {
        iterator tmp(*this);
        tmp -= i; return tmp;
    }

    difference_type operator-(iterator it) const
    {
        // for "x = it2 - it"
        assert(table == it.table);
        return pos - it.pos;
    }

    // Comparisons.
    bool operator==(const iterator& it) const
    {
        return table == it.table && pos == it.pos;
    }

    bool operator<(const iterator& it) const
    {
        assert(table == it.table);              // life is bad bad bad otherwise
        return pos < it.pos;
    }

    bool operator!=(const iterator& it) const { return !(*this == it); }
    bool operator<=(const iterator& it) const { return !(it < *this); }
    bool operator>(const iterator& it) const { return it < *this; }
    bool operator>=(const iterator& it) const { return !(*this < it); }

    // Here's the info we actually need to be an iterator
    tabletype *table;              // so we can dereference and bounds-check
    size_type pos;                 // index into the table
};

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
template <class tabletype>
class const_table_iterator
{
public:
    typedef table_iterator<tabletype> iterator;
    typedef const_table_iterator const_iterator;

    typedef std::random_access_iterator_tag iterator_category;
    typedef typename tabletype::value_type value_type;
    typedef typename tabletype::difference_type difference_type;
    typedef typename tabletype::size_type size_type;
    typedef typename tabletype::const_reference reference;  // we're const-only
    typedef typename tabletype::const_pointer pointer;

    // The "real" constructor
    const_table_iterator(const tabletype *tbl, size_type p)
        : table(tbl), pos(p) { }

    // The default constructor, used when I define vars of type table::iterator
    const_table_iterator() : table(NULL), pos(0) { }

    // The copy constructor, for when I say table::iterator foo = tbl.begin()
    // Also converts normal iterators to const iterators // not explicit on purpose
    const_table_iterator(const iterator &from)
        : table(from.table), pos(from.pos) { }

    // The default destructor is fine; we don't define one
    // The default operator= is fine; we don't define one

    // The main thing our iterator does is dereference.  If the table entry
    // we point to is empty, we return the default value type.
    reference operator*() const       { return (*table)[pos]; }
    pointer operator->() const        { return &(operator*()); }

    // Helper function to assert things are ok; eg pos is still in range
    void check() const
    {
        assert(table);
        assert(pos <= table->size());
    }

    // Arithmetic: we just do arithmetic on pos.  We don't even need to
    // do bounds checking, since STL doesn't consider that its job.  :-)
    const_iterator& operator+=(size_type t) { pos += t; check(); return *this; }
    const_iterator& operator-=(size_type t) { pos -= t; check(); return *this; }
    const_iterator& operator++()            { ++pos; check(); return *this; }
    const_iterator& operator--()            { --pos; check(); return *this; }
    const_iterator operator++(int)          
    {
        const_iterator tmp(*this); // for x++
        ++pos; check(); 
        return tmp; 
    }
    const_iterator operator--(int)          
    {
        const_iterator tmp(*this); // for x--
        --pos; check(); 
        return tmp;
    }
    const_iterator operator+(difference_type i) const
    {
        const_iterator tmp(*this);
        tmp += i;
        return tmp;
    }
    const_iterator operator-(difference_type i) const
    {
        const_iterator tmp(*this);
        tmp -= i;
        return tmp;
    }
    difference_type operator-(const_iterator it) const
    {
        // for "x = it2 - it"
        assert(table == it.table);
        return pos - it.pos;
    }
    reference operator[](difference_type n) const
    {
        return *(*this + n);            // simple though not totally efficient
    }

    // Comparisons.
    bool operator==(const const_iterator& it) const
    {
        return table == it.table && pos == it.pos;
    }

    bool operator<(const const_iterator& it) const
    {
        assert(table == it.table);              // life is bad bad bad otherwise
        return pos < it.pos;
    }
    bool operator!=(const const_iterator& it) const { return !(*this == it); }
    bool operator<=(const const_iterator& it) const { return !(it < *this); }
    bool operator>(const const_iterator& it) const { return it < *this; }
    bool operator>=(const const_iterator& it) const { return !(*this < it); }

    // Here's the info we actually need to be an iterator
    const tabletype *table;        // so we can dereference and bounds-check
    size_type pos;                 // index into the table
};

// ---------------------------------------------------------------------------
// This is a 2-D iterator.  You specify a begin and end over a list
// of *containers*.  We iterate over each container by iterating over
// it.  It's actually simple:
// VECTOR.begin() VECTOR[0].begin()  --------> VECTOR[0].end() ---,
//     |          ________________________________________________/
//     |          \_> VECTOR[1].begin()  -------->  VECTOR[1].end() -,
//     |          ___________________________________________________/
//     v          \_> ......
// VECTOR.end()
//
// It's impossible to do random access on one of these things in constant
// time, so it's just a bidirectional iterator.
//
// Unfortunately, because we need to use this for a non-empty iterator,
// we use ne_begin() and ne_end() instead of begin() and end()
// (though only going across, not down).
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
template <class T, class row_it, class col_it, class iter_type>
class Two_d_iterator
{
public:
    typedef Two_d_iterator iterator;
    typedef iter_type      iterator_category;
    typedef T              value_type;
    typedef std::ptrdiff_t difference_type;
    typedef T*             pointer;
    typedef T&             reference;

    explicit Two_d_iterator(row_it curr) : row_current(curr), col_current(0)
    {
        if (row_current && !row_current->is_marked())
        {
            col_current = row_current->ne_begin();
            advance_past_end();                 // in case cur->begin() == cur->end()
        }
    }

    explicit Two_d_iterator(row_it curr, col_it col) : row_current(curr), col_current(col)
    {
        assert(col);
    }

    // The default constructor
    Two_d_iterator() :  row_current(0), col_current(0) { }

    // Need this explicitly so we can convert normal iterators <=> const iterators
    // not explicit on purpose
    // ---------------------------------------------------------------------------
    template <class T2, class row_it2, class col_it2, class iter_type2>
    Two_d_iterator(const Two_d_iterator<T2, row_it2, col_it2, iter_type2>& it) :
        row_current (*(row_it *)&it.row_current),
        col_current (*(col_it *)&it.col_current)
    { }

    // The default destructor is fine; we don't define one
    // The default operator= is fine; we don't define one

    value_type& operator*() const  { return *(col_current); }
    value_type* operator->() const { return &(operator*()); }

    // Arithmetic: we just do arithmetic on pos.  We don't even need to
    // do bounds checking, since STL doesn't consider that its job.  :-)
    // NOTE: this is not amortized constant time!  What do we do about it?
    // ------------------------------------------------------------------
    void advance_past_end()
    {
        // used when col_current points to end()
        while (col_current == row_current->ne_end())
        {
            // end of current row
            // ------------------
            ++row_current;                                // go to beginning of next
            if (!row_current->is_marked())                // col is irrelevant at end
                col_current = row_current->ne_begin();
            else
                break;                                    // don't go past row_end
        }
    }

    friend size_t operator-(iterator l, iterator f)
    {
        if (f.row_current->is_marked())
            return 0;

        size_t diff(0);
        while (f != l)
        {
            ++diff;
            ++f;
        }
        return diff;
    }

    iterator& operator++()
    {
        // assert(!row_current->is_marked());               // how to ++ from there?
        ++col_current;
        advance_past_end();                              // in case col_current is at end()
        return *this;
    }

    iterator& operator--()
    {
        while (row_current->is_marked() ||
               col_current == row_current->ne_begin())
        {
            --row_current;
            col_current = row_current->ne_end();             // this is 1 too far
        }
        --col_current;
        return *this;
    }
    iterator operator++(int)       { iterator tmp(*this); ++*this; return tmp; }
    iterator operator--(int)       { iterator tmp(*this); --*this; return tmp; }


    // Comparisons.
    bool operator==(const iterator& it) const
    {
        return (row_current == it.row_current &&
                (!row_current || row_current->is_marked() || col_current == it.col_current));
    }

    bool operator!=(const iterator& it) const { return !(*this == it); }

    // Here's the info we actually need to be an iterator
    // These need to be public so we convert from iterator to const_iterator
    // ---------------------------------------------------------------------
    row_it row_current;
    col_it col_current;
};


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
template <class T, class row_it, class col_it, class iter_type, class Alloc>
class Two_d_destructive_iterator : public Two_d_iterator<T, row_it, col_it, iter_type>
{
public:
    typedef Two_d_destructive_iterator iterator;

    Two_d_destructive_iterator(Alloc &alloc, row_it curr) :
        _alloc(alloc)
    {
        this->row_current = curr;
        this->col_current = 0;
        if (this->row_current && !this->row_current->is_marked())
        {
            this->col_current = this->row_current->ne_begin();
            advance_past_end();                 // in case cur->begin() == cur->end()
        }
    }

    // Arithmetic: we just do arithmetic on pos.  We don't even need to
    // do bounds checking, since STL doesn't consider that its job.  :-)
    // NOTE: this is not amortized constant time!  What do we do about it?
    // ------------------------------------------------------------------
    void advance_past_end()
    {
        // used when col_current points to end()
        while (this->col_current == this->row_current->ne_end())
        {
            this->row_current->clear(_alloc, true);  // This is what differs from non-destructive iterators above

            // end of current row
            // ------------------
            ++this->row_current;                          // go to beginning of next
            if (!this->row_current->is_marked())          // col is irrelevant at end
                this->col_current = this->row_current->ne_begin();
            else
                break;                                    // don't go past row_end
        }
    }

    iterator& operator++()
    {
        // assert(!this->row_current->is_marked());         // how to ++ from there?
        ++this->col_current;
        advance_past_end();                              // in case col_current is at end()
        return *this;
    }

private:
    Two_d_destructive_iterator& operator=(const Two_d_destructive_iterator &o);

    Alloc &_alloc;
};


// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
#if defined(SPP_POPCNT_CHECK)
static inline bool spp_popcount_check()
{
    int cpuInfo[4] = { -1 };
    spp_cpuid(cpuInfo, 1);
    if (cpuInfo[2] & (1 << 23))
        return true;   // means SPP_POPCNT supported
    return false;
}
#endif

#if defined(SPP_POPCNT_CHECK) && defined(SPP_POPCNT)

static inline uint32_t spp_popcount(uint32_t i)
{
    static const bool s_ok = spp_popcount_check();
    return s_ok ? SPP_POPCNT(i) : s_spp_popcount_default(i);
}

#else

static inline uint32_t spp_popcount(uint32_t i)
{
#if defined(SPP_POPCNT)
    return static_cast<uint32_t>(SPP_POPCNT(i));
#else
    return s_spp_popcount_default(i);
#endif
}

#endif

#if defined(SPP_POPCNT_CHECK) && defined(SPP_POPCNT64)

static inline uint32_t spp_popcount(uint64_t i)
{
    static const bool s_ok = spp_popcount_check();
    return s_ok ? (uint32_t)SPP_POPCNT64(i) : s_spp_popcount_default(i);
}

#else

static inline uint32_t spp_popcount(uint64_t i)
{
#if defined(SPP_POPCNT64)
    return static_cast<uint32_t>(SPP_POPCNT64(i));
#elif 1
    return s_spp_popcount_default(i);
#endif
}

#endif

// ---------------------------------------------------------------------------
// SPARSE-TABLE
// ------------
// The idea is that a table with (logically) t buckets is divided
// into t/M *groups* of M buckets each.  (M is a constant, typically
// 32)  Each group is stored sparsely.
// Thus, inserting into the table causes some array to grow, which is
// slow but still constant time.  Lookup involves doing a
// logical-position-to-sparse-position lookup, which is also slow but
// constant time.  The larger M is, the slower these operations are
// but the less overhead (slightly).
//
// To store the sparse array, we store a bitmap B, where B[i] = 1 iff
// bucket i is non-empty.  Then to look up bucket i we really look up
// array[# of 1s before i in B].  This is constant time for fixed M.
//
// Terminology: the position of an item in the overall table (from
// 1 .. t) is called its "location."  The logical position in a group
// (from 1 .. M) is called its "position."  The actual location in
// the array (from 1 .. # of non-empty buckets in the group) is
// called its "offset."
// ---------------------------------------------------------------------------

template <class T, class Alloc>
class sparsegroup
{
public:
    // Basic types
    typedef T                                              value_type;
    typedef Alloc                                          allocator_type;
    typedef value_type&                                    reference;
    typedef const value_type&                              const_reference;
    typedef value_type*                                    pointer;
    typedef const value_type*                              const_pointer;

    typedef uint8_t                                        size_type;        // max # of buckets

    // These are our special iterators, that go over non-empty buckets in a
    // group.  These aren't const-only because you can change non-empty bcks.
    // ---------------------------------------------------------------------
    typedef pointer                                        ne_iterator;
    typedef const_pointer                                  const_ne_iterator;
    typedef std::reverse_iterator<ne_iterator>             reverse_ne_iterator;
    typedef std::reverse_iterator<const_ne_iterator>       const_reverse_ne_iterator;

    // We'll have versions for our special non-empty iterator too
    // ----------------------------------------------------------
    ne_iterator               ne_begin()         { return reinterpret_cast<pointer>(_group); }
    const_ne_iterator         ne_begin() const   { return reinterpret_cast<pointer>(_group); }
    const_ne_iterator         ne_cbegin() const  { return reinterpret_cast<pointer>(_group); }
    ne_iterator               ne_end()           { return reinterpret_cast<pointer>(_group + _num_items()); }
    const_ne_iterator         ne_end() const     { return reinterpret_cast<pointer>(_group + _num_items()); }
    const_ne_iterator         ne_cend() const    { return reinterpret_cast<pointer>(_group + _num_items()); }
    reverse_ne_iterator       ne_rbegin()        { return reverse_ne_iterator(ne_end()); }
    const_reverse_ne_iterator ne_rbegin() const  { return const_reverse_ne_iterator(ne_cend());  }
    const_reverse_ne_iterator ne_crbegin() const { return const_reverse_ne_iterator(ne_cend());  }
    reverse_ne_iterator       ne_rend()          { return reverse_ne_iterator(ne_begin()); }
    const_reverse_ne_iterator ne_rend() const    { return const_reverse_ne_iterator(ne_cbegin());  }
    const_reverse_ne_iterator ne_crend() const   { return const_reverse_ne_iterator(ne_cbegin());  }

private:
    // T can be std::pair<const K, V>, but sometime we need to cast to a mutable type
    // ------------------------------------------------------------------------------
    typedef typename spp_::cvt<T>::type                    mutable_value_type;
    typedef mutable_value_type &                           mutable_reference;
    typedef mutable_value_type *                           mutable_pointer;
    typedef const mutable_value_type *                     const_mutable_pointer;

    bool _bmtest(size_type i) const   { return !!(_bitmap & (static_cast<group_bm_type>(1) << i)); }
    void _bmset(size_type i)          { _bitmap |= static_cast<group_bm_type>(1) << i; }
    void _bmclear(size_type i)        { _bitmap &= ~(static_cast<group_bm_type>(1) << i); }

    bool _bme_test(size_type i) const { return !!(_bm_erased & (static_cast<group_bm_type>(1) << i)); }
    void _bme_set(size_type i)        { _bm_erased |= static_cast<group_bm_type>(1) << i; }
    void _bme_clear(size_type i)      { _bm_erased &= ~(static_cast<group_bm_type>(1) << i); }

    bool _bmtest_strict(size_type i) const
    { return !!((_bitmap | _bm_erased) & (static_cast<group_bm_type>(1) << i)); }


    static uint32_t _sizing(uint32_t n)
    {
#if !defined(SPP_ALLOC_SZ) || (SPP_ALLOC_SZ == 0)
        // aggressive allocation first, then decreasing as sparsegroups fill up
        // --------------------------------------------------------------------
        struct alloc_batch_size
        {
            // 32 bit bitmap
            // ........ .... .... .. .. .. .. .  .  .  .  .  .  .  .
            //     8     12   16  18 20 22 24 25 26   ...          32
            // ------------------------------------------------------
            SPP_CXX14_CONSTEXPR alloc_batch_size()
                : data()
            {
                uint8_t group_sz          = SPP_GROUP_SIZE / 4;
                uint8_t group_start_alloc = SPP_GROUP_SIZE / 8; //4;
                uint8_t alloc_sz          = group_start_alloc;
                for (int i=0; i<4; ++i)
                {
                    for (int j=0; j<group_sz; ++j)
                    {
                        if (j && j % group_start_alloc == 0)
                            alloc_sz += group_start_alloc;
                        data[i * group_sz + j] = alloc_sz;
                    }
                    if (group_start_alloc > 2)
                        group_start_alloc /= 2;
                    alloc_sz += group_start_alloc;
                }
            }
            uint8_t data[SPP_GROUP_SIZE];
        };

        static alloc_batch_size s_alloc_batch_sz;
        return n ? static_cast<uint32_t>(s_alloc_batch_sz.data[n-1]) : 0; // more aggressive alloc at the beginning

#elif (SPP_ALLOC_SZ == 1)
        // use as little memory as possible - slowest insert/delete in table
        // -----------------------------------------------------------------
        return n;
#else
        // decent compromise when SPP_ALLOC_SZ == 2
        // ----------------------------------------
        static size_type sz_minus_1 = SPP_ALLOC_SZ - 1;
        return (n + sz_minus_1) & ~sz_minus_1;
#endif
    }

    pointer _allocate_group(allocator_type &alloc, uint32_t n /* , bool tight = false */)
    {
        // ignore tight since we don't store num_alloc
        // num_alloc = (uint8_t)(tight ? n : _sizing(n));

        uint32_t num_alloc = (uint8_t)_sizing(n);
        _set_num_alloc(num_alloc);
        pointer retval = alloc.allocate(static_cast<size_type>(num_alloc));
        if (retval == NULL)
        {
            // the allocator is supposed to throw an exception if the allocation fails.
            throw_exception(std::bad_alloc());
        }
        return retval;
    }

    void _free_group(allocator_type &alloc, uint32_t num_alloc)
    {
        if (_group)
        {
            uint32_t num_buckets = _num_items();
            if (num_buckets)
            {
                mutable_pointer end_it = (mutable_pointer)(_group + num_buckets);
                for (mutable_pointer p = (mutable_pointer)_group; p != end_it; ++p)
                    p->~mutable_value_type();
            }
            alloc.deallocate(_group, (typename allocator_type::size_type)num_alloc);
            _group = NULL;
        }
    }

    // private because should not be called - no allocator!
    sparsegroup &operator=(const sparsegroup& x);

    static size_type _pos_to_offset(group_bm_type bm, size_type pos)
    {
        //return (size_type)((uint32_t)~((int32_t(-1) + pos) >> 31) & spp_popcount(bm << (SPP_GROUP_SIZE - pos)));
        //return (size_type)(pos ? spp_popcount(bm << (SPP_GROUP_SIZE - pos)) : 0);
        return static_cast<size_type>(spp_popcount(bm & ((static_cast<group_bm_type>(1) << pos) - 1)));
    }

public:

    // get_iter() in sparsetable needs it
    size_type pos_to_offset(size_type pos) const
    {
        return _pos_to_offset(_bitmap, pos);
    }

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4146)
#endif

    // Returns the (logical) position in the bm[] array, i, such that
    // bm[i] is the offset-th set bit in the array.  It is the inverse
    // of pos_to_offset.  get_pos() uses this function to find the index
    // of an ne_iterator in the table.  Bit-twiddling from
    // http://hackersdelight.org/basics.pdf
    // -----------------------------------------------------------------
    static size_type offset_to_pos(group_bm_type bm, size_type offset)
    {
        for (; offset > 0; offset--)
            bm &= (bm-1);  // remove right-most set bit

        // Clear all bits to the left of the rightmost bit (the &),
        // and then clear the rightmost bit but set all bits to the
        // right of it (the -1).
        // --------------------------------------------------------
        bm = (bm & -bm) - 1;
        return  static_cast<size_type>(spp_popcount(bm));
    }

#ifdef _MSC_VER
#pragma warning(pop)
#endif

    size_type offset_to_pos(size_type offset) const
    {
        return offset_to_pos(_bitmap, offset);
    }

public:
    // Constructors -- default and copy -- and destructor
    explicit sparsegroup() :
        _group(0), _bitmap(0), _bm_erased(0)
    {
        _set_num_items(0);
        _set_num_alloc(0);
    }

    sparsegroup(const sparsegroup& x) :
        _group(0), _bitmap(x._bitmap), _bm_erased(x._bm_erased)
    {
        _set_num_items(0);
        _set_num_alloc(0);
         assert(_group == 0); 
    }

    sparsegroup(const sparsegroup& x, allocator_type& a) :
        _group(0), _bitmap(x._bitmap), _bm_erased(x._bm_erased)
    {
        _set_num_items(0);
        _set_num_alloc(0);

        uint32_t num_items = x._num_items();
        if (num_items)
        {
            _group = _allocate_group(a, num_items /* , true */);
            _set_num_items(num_items);
            std::uninitialized_copy(x._group, x._group + num_items, _group);
        }
    }

    ~sparsegroup() { assert(_group == 0); }

    void destruct(allocator_type& a) { _free_group(a, _num_alloc()); }

    // Many STL algorithms use swap instead of copy constructors
    void swap(sparsegroup& x)
    {
        using std::swap;

        swap(_group, x._group);
        swap(_bitmap, x._bitmap);
        swap(_bm_erased, x._bm_erased);
#ifdef SPP_STORE_NUM_ITEMS
        swap(_num_buckets,   x._num_buckets);
        swap(_num_allocated, x._num_allocated);
#endif
    }

    // It's always nice to be able to clear a table without deallocating it
    void clear(allocator_type &alloc, bool erased)
    {
        _free_group(alloc, _num_alloc());
        _bitmap = 0;
        if (erased)
            _bm_erased = 0;
        _set_num_items(0);
        _set_num_alloc(0);
    }

    // Functions that tell you about size.  Alas, these aren't so useful
    // because our table is always fixed size.
    size_type size() const           { return static_cast<size_type>(SPP_GROUP_SIZE); }
    size_type max_size() const       { return static_cast<size_type>(SPP_GROUP_SIZE); }

    bool empty() const               { return false; }

    // We also may want to know how many *used* buckets there are
    size_type num_nonempty() const   { return (size_type)_num_items(); }

    // TODO(csilvers): make protected + friend
    // This is used by sparse_hashtable to get an element from the table
    // when we know it exists.
    reference unsafe_get(size_type i) const
    {
        // assert(_bmtest(i));
        return (reference)_group[pos_to_offset(i)];
    }

    typedef std::pair<pointer, bool> SetResult;

private:
    //typedef spp_::integral_constant<bool, spp_::is_relocatable<value_type>::value> check_relocatable;
    typedef spp_::true_type  realloc_ok_type;
    typedef spp_::false_type realloc_not_ok_type;

    //typedef spp_::zero_type  libc_reloc_type;
    //typedef spp_::one_type   spp_reloc_type;
    //typedef spp_::two_type   spp_not_reloc_type;
    //typedef spp_::three_type generic_alloc_type;

#if 1
    typedef typename if_<((spp_::is_same<allocator_type, libc_allocator<value_type> >::value ||
                           spp_::is_same<allocator_type,  spp_allocator<value_type> >::value) &&
                          spp_::is_relocatable<value_type>::value), realloc_ok_type, realloc_not_ok_type>::type
             check_alloc_type;
#else
    typedef typename if_<spp_::is_same<allocator_type, spp_allocator<value_type> >::value,
                         typename if_<spp_::is_relocatable<value_type>::value, spp_reloc_type, spp_not_reloc_type>::type,
                         typename if_<(spp_::is_same<allocator_type, libc_allocator<value_type> >::value &&
                                       spp_::is_relocatable<value_type>::value), libc_reloc_type, generic_alloc_type>::type >::type 
        check_alloc_type;
#endif


    //typedef if_<spp_::is_same<allocator_type, libc_allocator<value_type> >::value,
    //            libc_alloc_type,
    //            if_<spp_::is_same<allocator_type, spp_allocator<value_type> >::value,
    //                spp_alloc_type, user_alloc_type> > check_alloc_type;

    //typedef spp_::integral_constant<bool,
    //            (spp_::is_relocatable<value_type>::value &&
    //             (spp_::is_same<allocator_type, spp_allocator<value_type> >::value ||
    //              spp_::is_same<allocator_type, libc_allocator<value_type> >::value)) >
    //        realloc_and_memmove_ok;

    // ------------------------- memory at *p is uninitialized => need to construct
    void _init_val(mutable_value_type *p, reference val)
    {
#if !defined(SPP_NO_CXX11_RVALUE_REFERENCES)
        ::new (p) value_type(std::move((mutable_reference)val));
#else
        ::new (p) value_type((mutable_reference)val);
#endif
    }

    // ------------------------- memory at *p is uninitialized => need to construct
    void _init_val(mutable_value_type *p, const_reference val)
    {
        ::new (p) value_type(val);
    }

    // ------------------------------------------------ memory at *p is initialized
    void _set_val(value_type *p, reference val)
    {
#if !defined(SPP_NO_CXX11_RVALUE_REFERENCES)
        *(mutable_pointer)p = std::move((mutable_reference)val);
#else
        using std::swap;
        swap(*(mutable_pointer)p, *(mutable_pointer)&val);
#endif
    }

    // ------------------------------------------------ memory at *p is initialized
    void _set_val(value_type *p, const_reference val)
    {
        *(mutable_pointer)p = *(const_mutable_pointer)&val;
    }

    // Create space at _group[offset], assuming value_type is relocatable, and the 
    // allocator_type is the spp allocator.
    // return true if the slot was constructed (i.e. contains a valid value_type
    // ---------------------------------------------------------------------------------
    template <class Val>
    void _set_aux(allocator_type &alloc, size_type offset, Val &val, realloc_ok_type)
    {
        //static int x=0;  if (++x < 10) printf("x\n"); // check we are getting here

        uint32_t  num_items = _num_items();
        uint32_t  num_alloc = _sizing(num_items);

        if (num_items == num_alloc)
        {
            num_alloc = _sizing(num_items + 1);
            _group = alloc.reallocate(_group, num_alloc);
            _set_num_alloc(num_alloc);
        }

        for (uint32_t i = num_items; i > offset; --i)
            memcpy(static_cast<void *>(_group + i), _group + i-1, sizeof(*_group));

        _init_val((mutable_pointer)(_group + offset), val);
    }

    // Create space at _group[offset], assuming value_type is *not* relocatable, and the 
    // allocator_type is the spp allocator.
    // return true if the slot was constructed (i.e. contains a valid value_type
    // ---------------------------------------------------------------------------------
    template <class Val>
    void _set_aux(allocator_type &alloc, size_type offset, Val &val, realloc_not_ok_type)
    {
        uint32_t  num_items = _num_items();
        uint32_t  num_alloc = _sizing(num_items);

        //assert(num_alloc == (uint32_t)_num_allocated);
        if (num_items < num_alloc)
        {
            // create new object at end and rotate it to position
            _init_val((mutable_pointer)&_group[num_items], val);
            std::rotate((mutable_pointer)(_group + offset),
                        (mutable_pointer)(_group + num_items),
                        (mutable_pointer)(_group + num_items + 1));
            return;
        }

        // This is valid because 0 <= offset <= num_items
        pointer p = _allocate_group(alloc, _sizing(num_items + 1));
        if (offset)
            std::uninitialized_copy(MK_MOVE_IT((mutable_pointer)_group),
                                    MK_MOVE_IT((mutable_pointer)(_group + offset)),
                                    (mutable_pointer)p);
        if (num_items > offset)
            std::uninitialized_copy(MK_MOVE_IT((mutable_pointer)(_group + offset)),
                                    MK_MOVE_IT((mutable_pointer)(_group + num_items)),
                                    (mutable_pointer)(p + offset + 1));
        _init_val((mutable_pointer)(p + offset), val);
        _free_group(alloc, num_alloc);
        _group = p;
    }

    // ----------------------------------------------------------------------------------
    template <class Val>
    void _set(allocator_type &alloc, size_type i, size_type offset, Val &val)
    {
        if (!_bmtest(i))
        {
            _set_aux(alloc, offset, val, check_alloc_type());
            _incr_num_items();
            _bmset(i);
        }
        else
            _set_val(&_group[offset], val);
    }

public:

    // This returns the pointer to the inserted item
    // ---------------------------------------------
    template <class Val>
    pointer set(allocator_type &alloc, size_type i, Val &val)
    {
        _bme_clear(i); // in case this was an "erased" location

        size_type offset = pos_to_offset(i);
        _set(alloc, i, offset, val);            // may change _group pointer
        return (pointer)(_group + offset);
    }

    // We let you see if a bucket is non-empty without retrieving it
    // -------------------------------------------------------------
    bool test(size_type i) const
    {
        return _bmtest(i);
    }

    // also tests for erased values
    // ----------------------------
    bool test_strict(size_type i) const
    {
        return _bmtest_strict(i);
    }

private:
    // Shrink the array, assuming value_type is relocatable, and the 
    // allocator_type is the libc allocator (supporting reallocate).
    // -------------------------------------------------------------
    void _group_erase_aux(allocator_type &alloc, size_type offset, realloc_ok_type)
    {
        // static int x=0;  if (++x < 10) printf("Y\n"); // check we are getting here
        uint32_t  num_items = _num_items();
        uint32_t  num_alloc = _sizing(num_items);

        if (num_items == 1)
        {
            assert(offset == 0);
            _free_group(alloc, num_alloc);
            _set_num_alloc(0);
            return;
        }

        _group[offset].~value_type();

        for (size_type i = offset; i < num_items - 1; ++i)
            memcpy(static_cast<void *>(_group + i), _group + i + 1, sizeof(*_group));

        if (_sizing(num_items - 1) != num_alloc)
        {
            num_alloc = _sizing(num_items - 1);
            assert(num_alloc);            // because we have at least 1 item left
            _set_num_alloc(num_alloc);
            _group = alloc.reallocate(_group, num_alloc);
        }
    }

    // Shrink the array, without any special assumptions about value_type and
    // allocator_type.
    // --------------------------------------------------------------------------
    void _group_erase_aux(allocator_type &alloc, size_type offset, realloc_not_ok_type)
    {
        uint32_t  num_items = _num_items();
        uint32_t  num_alloc   = _sizing(num_items);

        if (_sizing(num_items - 1) != num_alloc)
        {
            pointer p = 0;
            if (num_items > 1)
            {
                p = _allocate_group(alloc, num_items - 1);
                if (offset)
                    std::uninitialized_copy(MK_MOVE_IT((mutable_pointer)(_group)),
                                            MK_MOVE_IT((mutable_pointer)(_group + offset)),
                                            (mutable_pointer)(p));
                if (static_cast<uint32_t>(offset + 1) < num_items)
                    std::uninitialized_copy(MK_MOVE_IT((mutable_pointer)(_group + offset + 1)),
                                            MK_MOVE_IT((mutable_pointer)(_group + num_items)),
                                            (mutable_pointer)(p + offset));
            }
            else
            {
                assert(offset == 0);
                _set_num_alloc(0);
            }
            _free_group(alloc, num_alloc);
            _group = p;
        }
        else
        {
            std::rotate((mutable_pointer)(_group + offset),
                        (mutable_pointer)(_group + offset + 1),
                        (mutable_pointer)(_group + num_items));
            ((mutable_pointer)(_group + num_items - 1))->~mutable_value_type();
        }
    }

    void _group_erase(allocator_type &alloc, size_type offset)
    {
        _group_erase_aux(alloc, offset, check_alloc_type());
    }

public:
    template <class twod_iter>
    bool erase_ne(allocator_type &alloc, twod_iter &it)
    {
        assert(_group && it.col_current != ne_end());
        size_type offset = (size_type)(it.col_current - ne_begin());
        size_type pos    = offset_to_pos(offset);

        if (_num_items() <= 1)
        {
            clear(alloc, false);
            it.col_current = 0;
        }
        else
        {
            _group_erase(alloc, offset);
            _decr_num_items();
            _bmclear(pos);

            // in case _group_erase reallocated the buffer
            it.col_current = reinterpret_cast<pointer>(_group) + offset;
        }
        _bme_set(pos);  // remember that this position has been erased
        it.advance_past_end();
        return true;
    }


    // This takes the specified elements out of the group.  This is
    // "undefining", rather than "clearing".
    // TODO(austern): Make this exception safe: handle exceptions from
    // value_type's copy constructor.
    // ---------------------------------------------------------------
    void erase(allocator_type &alloc, size_type i)
    {
        if (_bmtest(i))
        {
            // trivial to erase empty bucket
            if (_num_items() == 1)
                clear(alloc, false);
            else
            {
                _group_erase(alloc, pos_to_offset(i));
                _decr_num_items();
                _bmclear(i);
            }
            _bme_set(i); // remember that this position has been erased
        }
    }

    // I/O
    // We support reading and writing groups to disk.  We don't store
    // the actual array contents (which we don't know how to store),
    // just the bitmap and size.  Meant to be used with table I/O.
    // --------------------------------------------------------------
    template <typename OUTPUT> bool write_metadata(OUTPUT *fp) const
    {
        // warning: we write 4 or 8 bytes for the bitmap, instead of 6 in the
        //          original google sparsehash
        // ------------------------------------------------------------------
        if (!sparsehash_internal::write_data(fp, &_bitmap, sizeof(_bitmap)))
            return false;

        return true;
    }

    // Reading destroys the old group contents!  Returns true if all was ok.
    template <typename INPUT> bool read_metadata(allocator_type &alloc, INPUT *fp)
    {
        clear(alloc, true);

        if (!sparsehash_internal::read_data(fp, &_bitmap, sizeof(_bitmap)))
            return false;

        // We'll allocate the space, but we won't fill it: it will be
        // left as uninitialized raw memory.
        uint32_t num_items = spp_popcount(_bitmap); // yes, _num_buckets not set
        _set_num_items(num_items);
        _group = num_items ? _allocate_group(alloc, num_items/* , true */) : 0;
        return true;
    }

    // Again, only meaningful if value_type is a POD.
    template <typename INPUT> bool read_nopointer_data(INPUT *fp)
    {
        for (ne_iterator it = ne_begin(); it != ne_end(); ++it)
            if (!sparsehash_internal::read_data(fp, &(*it), sizeof(*it)))
                return false;
        return true;
    }

    // If your keys and values are simple enough, we can write them
    // to disk for you.  "simple enough" means POD and no pointers.
    // However, we don't try to normalize endianness.
    // ------------------------------------------------------------
    template <typename OUTPUT> bool write_nopointer_data(OUTPUT *fp) const
    {
        for (const_ne_iterator it = ne_begin(); it != ne_end(); ++it)
            if (!sparsehash_internal::write_data(fp, &(*it), sizeof(*it)))
                return false;
        return true;
    }


    // Comparisons.  We only need to define == and < -- we get
    // != > <= >= via relops.h (which we happily included above).
    // Note the comparisons are pretty arbitrary: we compare
    // values of the first index that isn't equal (using default
    // value for empty buckets).
    // ---------------------------------------------------------
    bool operator==(const sparsegroup& x) const
    {
        return (_bitmap == x._bitmap &&
                _bm_erased == x._bm_erased &&
                std::equal(_group, _group + _num_items(), x._group));
    }

    bool operator<(const sparsegroup& x) const
    {
        // also from <algorithm>
        return std::lexicographical_compare(_group, _group + _num_items(),
                                            x._group, x._group + x._num_items());
    }

    bool operator!=(const sparsegroup& x) const { return !(*this == x); }
    bool operator<=(const sparsegroup& x) const { return !(x < *this); }
    bool operator> (const sparsegroup& x) const { return x < *this; }
    bool operator>=(const sparsegroup& x) const { return !(*this < x); }

    void mark()            { _group = (value_type *)static_cast<uintptr_t>(-1); }
    bool is_marked() const { return _group == (value_type *)static_cast<uintptr_t>(-1); }

private:
    // ---------------------------------------------------------------------------
    template <class A>
    class alloc_impl : public A
    {
    public:
        typedef typename A::pointer pointer;
        typedef typename A::size_type size_type;

        // Convert a normal allocator to one that has realloc_or_die()
        explicit alloc_impl(const A& a) : A(a) { }

        // realloc_or_die should only be used when using the default
        // allocator (spp::spp_allocator).
        pointer realloc_or_die(pointer /*ptr*/, size_type /*n*/)
        {
            throw_exception(std::runtime_error("realloc_or_die is only supported for spp::spp_allocator\n"));
            return NULL;
        }
    };

    // A template specialization of alloc_impl for
    // spp::libc_allocator that can handle realloc_or_die.
    // -----------------------------------------------------------
    template <class A>
    class alloc_impl<spp_::libc_allocator<A> > : public spp_::libc_allocator<A>
    {
    public:
        typedef typename spp_::libc_allocator<A>::pointer pointer;
        typedef typename spp_::libc_allocator<A>::size_type size_type;

        explicit alloc_impl(const spp_::libc_allocator<A>& a)
            : spp_::libc_allocator<A>(a)
        { }

        pointer realloc_or_die(pointer ptr, size_type n)
        {
            pointer retval = this->reallocate(ptr, n);
            if (retval == NULL) 
            {
                // the allocator is supposed to throw an exception if the allocation fails.
                throw_exception(std::bad_alloc());
            }
            return retval;
        }
    };

    // A template specialization of alloc_impl for
    // spp::spp_allocator that can handle realloc_or_die.
    // -----------------------------------------------------------
    template <class A>
    class alloc_impl<spp_::spp_allocator<A> > : public spp_::spp_allocator<A>
    {
    public:
        typedef typename spp_::spp_allocator<A>::pointer pointer;
        typedef typename spp_::spp_allocator<A>::size_type size_type;

        explicit alloc_impl(const spp_::spp_allocator<A>& a)
            : spp_::spp_allocator<A>(a)
        { }

        pointer realloc_or_die(pointer ptr, size_type n)
        {
            pointer retval = this->reallocate(ptr, n);
            if (retval == NULL) 
            {
                // the allocator is supposed to throw an exception if the allocation fails.
                throw_exception(std::bad_alloc());
            }
            return retval;
        }
    };


#ifdef SPP_STORE_NUM_ITEMS
    uint32_t _num_items() const           { return (uint32_t)_num_buckets; }
    void     _set_num_items(uint32_t val) { _num_buckets = static_cast<size_type>(val); }
    void     _incr_num_items()            { ++_num_buckets; }
    void     _decr_num_items()            { --_num_buckets; }
    uint32_t _num_alloc() const           { return (uint32_t)_num_allocated; }
    void     _set_num_alloc(uint32_t val) { _num_allocated = static_cast<size_type>(val); }
#else
    uint32_t _num_items() const           { return spp_popcount(_bitmap); }
    void     _set_num_items(uint32_t )    { }
    void     _incr_num_items()            { }
    void     _decr_num_items()            { }
    uint32_t _num_alloc() const           { return _sizing(_num_items()); }
    void     _set_num_alloc(uint32_t val) { }
#endif

    // The actual data
    // ---------------
    value_type *         _group;                             // (small) array of T's
    group_bm_type        _bitmap;
    group_bm_type        _bm_erased;                         // ones where items have been erased

#ifdef SPP_STORE_NUM_ITEMS
    size_type            _num_buckets;
    size_type            _num_allocated;
#endif
};

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
template <class T, class Alloc>
class sparsetable
{
public:
    typedef T                                             value_type;
    typedef Alloc                                         allocator_type;
    typedef sparsegroup<value_type, allocator_type>       group_type;

private:
    typedef typename Alloc::template rebind<group_type>::other group_alloc_type;
    typedef typename group_alloc_type::size_type          group_size_type;

public:
    // Basic types
    // -----------
    typedef typename allocator_type::size_type            size_type;
    typedef typename allocator_type::difference_type      difference_type;
    typedef value_type&                                   reference;
    typedef const value_type&                             const_reference;
    typedef value_type*                                   pointer;
    typedef const value_type*                             const_pointer;

    typedef group_type&                                   GroupsReference;
    typedef const group_type&                             GroupsConstReference;

    typedef typename group_type::ne_iterator              ColIterator;
    typedef typename group_type::const_ne_iterator        ColConstIterator;

    typedef table_iterator<sparsetable<T, allocator_type> >        iterator;       // defined with index
    typedef const_table_iterator<sparsetable<T, allocator_type> >  const_iterator; // defined with index
    typedef std::reverse_iterator<const_iterator>         const_reverse_iterator;
    typedef std::reverse_iterator<iterator>               reverse_iterator;

    // These are our special iterators, that go over non-empty buckets in a
    // table.  These aren't const only because you can change non-empty bcks.
    // ----------------------------------------------------------------------
    typedef Two_d_iterator<T,
                           group_type *,
                           ColIterator,
                           std::bidirectional_iterator_tag> ne_iterator;

    typedef Two_d_iterator<const T,
                           const group_type *,
                           ColConstIterator,
                           std::bidirectional_iterator_tag> const_ne_iterator;

    // Another special iterator: it frees memory as it iterates (used to resize).
    // Obviously, you can only iterate over it once, which is why it's an input iterator
    // ---------------------------------------------------------------------------------
    typedef Two_d_destructive_iterator<T,
                                       group_type *,
                                       ColIterator,
                                       std::input_iterator_tag,
                                       allocator_type>     destructive_iterator;

    typedef std::reverse_iterator<ne_iterator>               reverse_ne_iterator;
    typedef std::reverse_iterator<const_ne_iterator>         const_reverse_ne_iterator;


    // Iterator functions
    // ------------------
    iterator               begin()         { return iterator(this, 0); }
    const_iterator         begin() const   { return const_iterator(this, 0); }
    const_iterator         cbegin() const  { return const_iterator(this, 0); }
    iterator               end()           { return iterator(this, size()); }
    const_iterator         end() const     { return const_iterator(this, size()); }
    const_iterator         cend() const    { return const_iterator(this, size()); }
    reverse_iterator       rbegin()        { return reverse_iterator(end()); }
    const_reverse_iterator rbegin() const  { return const_reverse_iterator(cend()); }
    const_reverse_iterator crbegin() const { return const_reverse_iterator(cend()); }
    reverse_iterator       rend()          { return reverse_iterator(begin()); }
    const_reverse_iterator rend() const    { return const_reverse_iterator(cbegin()); }
    const_reverse_iterator crend() const   { return const_reverse_iterator(cbegin()); }

    // Versions for our special non-empty iterator
    // ------------------------------------------
    ne_iterator       ne_begin()           { return ne_iterator      (_first_group); }
    const_ne_iterator ne_begin() const     { return const_ne_iterator(_first_group); }
    const_ne_iterator ne_cbegin() const    { return const_ne_iterator(_first_group); }
    ne_iterator       ne_end()             { return ne_iterator      (_last_group); }
    const_ne_iterator ne_end() const       { return const_ne_iterator(_last_group); }
    const_ne_iterator ne_cend() const      { return const_ne_iterator(_last_group); }

    reverse_ne_iterator       ne_rbegin()        { return reverse_ne_iterator(ne_end()); }
    const_reverse_ne_iterator ne_rbegin() const  { return const_reverse_ne_iterator(ne_end());  }
    const_reverse_ne_iterator ne_crbegin() const { return const_reverse_ne_iterator(ne_end());  }
    reverse_ne_iterator       ne_rend()          { return reverse_ne_iterator(ne_begin()); }
    const_reverse_ne_iterator ne_rend() const    { return const_reverse_ne_iterator(ne_begin()); }
    const_reverse_ne_iterator ne_crend() const   { return const_reverse_ne_iterator(ne_begin()); }

    destructive_iterator destructive_begin()
    {
        return destructive_iterator(_alloc, _first_group);
    }

    destructive_iterator destructive_end()
    {
        return destructive_iterator(_alloc, _last_group);
    }

    // How to deal with the proper group
    static group_size_type num_groups(size_type num)
    {
        // how many to hold num buckets
        return num == 0 ? (group_size_type)0 :
            (group_size_type)(((num-1) / SPP_GROUP_SIZE) + 1);
    }

    typename group_type::size_type pos_in_group(size_type i) const
    {
        return static_cast<typename group_type::size_type>(i & SPP_MASK_);
    }

    size_type group_num(size_type i) const
    {
        return (size_type)(i >> SPP_SHIFT_);
    }

    GroupsReference which_group(size_type i)
    {
        return _first_group[group_num(i)];
    }

    GroupsConstReference which_group(size_type i) const
    {
        return _first_group[group_num(i)];
    }

    void _alloc_group_array(group_size_type sz, group_type *&first, group_type *&last)
    {
        if (sz)
        {
            first = _group_alloc.allocate((size_type)(sz + 1)); // + 1 for end marker
            first[sz].mark();                      // for the ne_iterator
            last = first + sz;
        }
    }

    void _free_group_array(group_type *&first, group_type *&last)
    {
        if (first)
        {
            _group_alloc.deallocate(first, (group_size_type)(last - first + 1)); // + 1 for end marker
            first = last = 0;
        }
    }

    void _allocate_groups(size_type sz)
    {
        if (sz)
        {
            _alloc_group_array(sz, _first_group, _last_group);
            std::uninitialized_fill(_first_group, _last_group, group_type());
        }
    }

    void _free_groups()
    {
        if (_first_group)
        {
            for (group_type *g = _first_group; g != _last_group; ++g)
                g->destruct(_alloc);
            _free_group_array(_first_group, _last_group);
        }
    }

    void _cleanup()
    {
        _free_groups();    // sets _first_group = _last_group = 0
        _table_size  = 0;
        _num_buckets = 0;
    }

    void _init()
    {
        _first_group = 0;
        _last_group  = 0;
        _table_size  = 0;
        _num_buckets = 0;
    }

    void _copy(const sparsetable &o)
    {
        _table_size = o._table_size;
        _num_buckets = o._num_buckets;
        _alloc = o._alloc;                // todo - copy or move allocator according to...
        _group_alloc = o._group_alloc;    // http://en.cppreference.com/w/cpp/container/unordered_map/unordered_map

        group_size_type sz = (group_size_type)(o._last_group - o._first_group);
        if (sz)
        {
            _alloc_group_array(sz, _first_group, _last_group);
            for (group_size_type i=0; i<sz; ++i)
                new (_first_group + i) group_type(o._first_group[i], _alloc);
        }
    }

public:
    // Constructors -- default, normal (when you specify size), and copy
    explicit sparsetable(size_type sz = 0, const allocator_type &alloc = allocator_type()) :
        _first_group(0),
        _last_group(0),
        _table_size(sz),
        _num_buckets(0),
        _group_alloc(alloc),
        _alloc(alloc)
                       // todo - copy or move allocator according to
                       // http://en.cppreference.com/w/cpp/container/unordered_map/unordered_map
    {
        _allocate_groups(num_groups(sz));
    }

    ~sparsetable()
    {
        _free_groups();
    }

    sparsetable(const sparsetable &o)
    {
        _init();
        _copy(o);
    }

    sparsetable& operator=(const sparsetable &o)
    {
        _cleanup();
        _copy(o);
        return *this;
    }


#if !defined(SPP_NO_CXX11_RVALUE_REFERENCES)
    sparsetable(sparsetable&& o)
    {
        _init();
        this->swap(o);
    }

    sparsetable(sparsetable&& o, const allocator_type &alloc)
    {
        _init();
        this->swap(o);
        _alloc = alloc; // [gp todo] is this correct?
    }

    sparsetable& operator=(sparsetable&& o)
    {
        _cleanup();
        this->swap(o);
        return *this;
    }
#endif

    // Many STL algorithms use swap instead of copy constructors
    void swap(sparsetable& o)
    {
        using std::swap;

        swap(_first_group, o._first_group);
        swap(_last_group,  o._last_group);
        swap(_table_size,  o._table_size);
        swap(_num_buckets, o._num_buckets);
        if (_alloc != o._alloc)
            swap(_alloc, o._alloc);
        if (_group_alloc != o._group_alloc)
            swap(_group_alloc, o._group_alloc);
    }

    // It's always nice to be able to clear a table without deallocating it
    void clear()
    {
        _free_groups();
        _num_buckets = 0;
        _table_size = 0;
    }

    inline allocator_type get_allocator() const
    {
        return _alloc;
    }


    // Functions that tell you about size.
    // NOTE: empty() is non-intuitive!  It does not tell you the number
    // of not-empty buckets (use num_nonempty() for that).  Instead
    // it says whether you've allocated any buckets or not.
    // ----------------------------------------------------------------
    size_type size() const           { return _table_size; }
    size_type max_size() const       { return _alloc.max_size(); }
    bool empty() const               { return _table_size == 0; }
    size_type num_nonempty() const   { return _num_buckets; }

    // OK, we'll let you resize one of these puppies
    void resize(size_type new_size)
    {
        group_size_type sz = num_groups(new_size);
        group_size_type old_sz = (group_size_type)(_last_group - _first_group);

        if (sz != old_sz)
        {
            // resize group array
            // ------------------
            group_type *first = 0, *last = 0;
            if (sz)
            {
                _alloc_group_array(sz, first, last);
                memcpy(static_cast<void *>(first), _first_group, sizeof(*first) * (std::min)(sz, old_sz));
            }

            if (sz < old_sz)
            {
                for (group_type *g = _first_group + sz; g != _last_group; ++g)
                    g->destruct(_alloc);
            }
            else
                std::uninitialized_fill(first + old_sz, last, group_type());

            _free_group_array(_first_group, _last_group);
            _first_group = first;
            _last_group  = last;
        }
#if 0
        // used only in test program
        // todo: fix if sparsetable to be used directly
        // --------------------------------------------
        if (new_size < _table_size)
        {
            // lower num_buckets, clear last group
            if (pos_in_group(new_size) > 0)     // need to clear inside last group
                groups.back().erase(_alloc, groups.back().begin() + pos_in_group(new_size),
                                    groups.back().end());
            _num_buckets = 0;                   // refigure # of used buckets
            for (const group_type *group = _first_group; group != _last_group; ++group)
                _num_buckets += group->num_nonempty();
        }
#endif
        _table_size = new_size;
    }

    // We let you see if a bucket is non-empty without retrieving it
    // -------------------------------------------------------------
    bool test(size_type i) const
    {
        // assert(i < _table_size);
        return which_group(i).test(pos_in_group(i));
    }

    // also tests for erased values
    // ----------------------------
    bool test_strict(size_type i) const
    {
        // assert(i < _table_size);
        return which_group(i).test_strict(pos_in_group(i));
    }

    friend struct GrpPos;

    struct GrpPos
    {
        typedef typename sparsetable::ne_iterator ne_iter;
        GrpPos(const sparsetable &table, size_type i) :
            grp(table.which_group(i)), pos(table.pos_in_group(i)) {}

        bool test_strict() const { return grp.test_strict(pos); }
        bool test() const        { return grp.test(pos); }
        typename sparsetable::reference unsafe_get() const { return  grp.unsafe_get(pos); }
        ne_iter get_iter(typename sparsetable::reference ref)
        {
            return ne_iter((group_type *)&grp, &ref);
        }

        void erase(sparsetable &table) // item *must* be present
        {
            assert(table._num_buckets);
            ((group_type &)grp).erase(table._alloc, pos);
            --table._num_buckets;
        }

    private:
        GrpPos* operator=(const GrpPos&);

        const group_type &grp;
        typename group_type::size_type pos;
    };

    bool test(iterator pos) const
    {
        return which_group(pos.pos).test(pos_in_group(pos.pos));
    }

    bool test(const_iterator pos) const
    {
        return which_group(pos.pos).test(pos_in_group(pos.pos));
    }

    // TODO(csilvers): make protected + friend
    // This is used by sparse_hashtable to get an element from the table
    // when we know it exists (because the caller has called test(i)).
    // -----------------------------------------------------------------
    reference unsafe_get(size_type i) const
    {
        assert(i < _table_size);
        // assert(test(i));
        return which_group(i).unsafe_get(pos_in_group(i));
    }

    // Needed for hashtables, gets as a ne_iterator.  Crashes for empty bcks
    const_ne_iterator get_iter(size_type i) const
    {
        //assert(test(i));    // how can a ne_iterator point to an empty bucket?

        size_type grp_idx = group_num(i);

        return const_ne_iterator(_first_group + grp_idx,
                                 (_first_group[grp_idx].ne_begin() +
                                  _first_group[grp_idx].pos_to_offset(pos_in_group(i))));
    }

    const_ne_iterator get_iter(size_type i, ColIterator col_it) const
    {
        return const_ne_iterator(_first_group + group_num(i), col_it);
    }

    // For nonempty we can return a non-const version
    ne_iterator get_iter(size_type i)
    {
        //assert(test(i));    // how can a nonempty_iterator point to an empty bucket?

        size_type grp_idx = group_num(i);

        return ne_iterator(_first_group + grp_idx,
                           (_first_group[grp_idx].ne_begin() +
                            _first_group[grp_idx].pos_to_offset(pos_in_group(i))));
    }

    ne_iterator get_iter(size_type i, ColIterator col_it)
    {
        return ne_iterator(_first_group + group_num(i), col_it);
    }

    // And the reverse transformation.
    size_type get_pos(const const_ne_iterator& it) const
    {
        difference_type current_row = it.row_current - _first_group;
        difference_type current_col = (it.col_current - _first_group[current_row].ne_begin());
        return ((current_row * SPP_GROUP_SIZE) +
                _first_group[current_row].offset_to_pos(current_col));
    }

    // Val can be reference or const_reference
    // ---------------------------------------
    template <class Val>
    reference set(size_type i, Val &val)
    {
        assert(i < _table_size);
        group_type &group = which_group(i);
        typename group_type::size_type old_numbuckets = group.num_nonempty();
        pointer p(group.set(_alloc, pos_in_group(i), val));
        _num_buckets += group.num_nonempty() - old_numbuckets;
        return *p;
    }

    // used in _move_from (where we can move the old value instead of copying it
    void move(size_type i, reference val)
    {
        assert(i < _table_size);
        which_group(i).set(_alloc, pos_in_group(i), val);
        ++_num_buckets;
    }

    // This takes the specified elements out of the table.
    // --------------------------------------------------
    void erase(size_type i)
    {
        assert(i < _table_size);

        GroupsReference grp(which_group(i));
        typename group_type::size_type old_numbuckets = grp.num_nonempty();
        grp.erase(_alloc, pos_in_group(i));
        _num_buckets += grp.num_nonempty() - old_numbuckets;
    }

    void erase(iterator pos)
    {
        erase(pos.pos);
    }

    void erase(iterator start_it, iterator end_it)
    {
        // This could be more efficient, but then we'd need to figure
        // out if we spanned groups or not.  Doesn't seem worth it.
        for (; start_it != end_it; ++start_it)
            erase(start_it);
    }

    const_ne_iterator erase(const_ne_iterator it)
    {
        ne_iterator res(it);
        if (res.row_current->erase_ne(_alloc, res))
            _num_buckets--;
        return res;
    }

    const_ne_iterator erase(const_ne_iterator f, const_ne_iterator l)
    {
        size_t diff = l - f;
        while (diff--)
            f = erase(f);
        return f;
    }

    // We support reading and writing tables to disk.  We don't store
    // the actual array contents (which we don't know how to store),
    // just the groups and sizes.  Returns true if all went ok.

private:
    // Every time the disk format changes, this should probably change too
    typedef unsigned long MagicNumberType;
    static const MagicNumberType MAGIC_NUMBER = 0x24687531;

    // Old versions of this code write all data in 32 bits.  We need to
    // support these files as well as having support for 64-bit systems.
    // So we use the following encoding scheme: for values < 2^32-1, we
    // store in 4 bytes in big-endian order.  For values > 2^32, we
    // store 0xFFFFFFF followed by 8 bytes in big-endian order.  This
    // causes us to mis-read old-version code that stores exactly
    // 0xFFFFFFF, but I don't think that is likely to have happened for
    // these particular values.
    template <typename OUTPUT, typename IntType>
    static bool write_32_or_64(OUTPUT* fp, IntType value)
    {
        if (value < 0xFFFFFFFFULL)        // fits in 4 bytes
        {
            if (!sparsehash_internal::write_bigendian_number(fp, value, 4))
                return false;
        }
        else
        {
            if (!sparsehash_internal::write_bigendian_number(fp, 0xFFFFFFFFUL, 4))
                return false;
            if (!sparsehash_internal::write_bigendian_number(fp, value, 8))
                return false;
        }
        return true;
    }

    template <typename INPUT, typename IntType>
    static bool read_32_or_64(INPUT* fp, IntType *value)
    {
        // reads into value
        MagicNumberType first4 = 0;   // a convenient 32-bit unsigned type
        if (!sparsehash_internal::read_bigendian_number(fp, &first4, 4))
            return false;

        if (first4 < 0xFFFFFFFFULL)
        {
            *value = first4;
        }
        else
        {
            if (!sparsehash_internal::read_bigendian_number(fp, value, 8))
                return false;
        }
        return true;
    }

public:
    // read/write_metadata() and read_write/nopointer_data() are DEPRECATED.
    // Use serialize() and unserialize(), below, for new code.

    template <typename OUTPUT>
    bool write_metadata(OUTPUT *fp) const
    {
        if (!write_32_or_64(fp, MAGIC_NUMBER))  return false;
        if (!write_32_or_64(fp, _table_size))  return false;
        if (!write_32_or_64(fp, _num_buckets))  return false;

        for (const group_type *group = _first_group; group != _last_group; ++group)
            if (group->write_metadata(fp) == false)
                return false;
        return true;
    }

    // Reading destroys the old table contents!  Returns true if read ok.
    template <typename INPUT>
    bool read_metadata(INPUT *fp)
    {
        size_type magic_read = 0;
        if (!read_32_or_64(fp, &magic_read))  return false;
        if (magic_read != MAGIC_NUMBER)
        {
            clear();                        // just to be consistent
            return false;
        }

        if (!read_32_or_64(fp, &_table_size))  return false;
        if (!read_32_or_64(fp, &_num_buckets))  return false;

        resize(_table_size);                    // so the vector's sized ok
        for (group_type *group = _first_group; group != _last_group; ++group)
            if (group->read_metadata(_alloc, fp) == false)
                return false;
        return true;
    }

    // This code is identical to that for SparseGroup
    // If your keys and values are simple enough, we can write them
    // to disk for you.  "simple enough" means no pointers.
    // However, we don't try to normalize endianness
    bool write_nopointer_data(FILE *fp) const
    {
        for (const_ne_iterator it = ne_begin(); it != ne_end(); ++it)
            if (!fwrite(&*it, sizeof(*it), 1, fp))
                return false;
        return true;
    }

    // When reading, we have to override the potential const-ness of *it
    bool read_nopointer_data(FILE *fp)
    {
        for (ne_iterator it = ne_begin(); it != ne_end(); ++it)
            if (!fread(reinterpret_cast<void*>(&(*it)), sizeof(*it), 1, fp))
                return false;
        return true;
    }

    // INPUT and OUTPUT must be either a FILE, *or* a C++ stream
    //    (istream, ostream, etc) *or* a class providing
    //    Read(void*, size_t) and Write(const void*, size_t)
    //    (respectively), which writes a buffer into a stream
    //    (which the INPUT/OUTPUT instance presumably owns).

    typedef sparsehash_internal::pod_serializer<value_type> NopointerSerializer;

    // ValueSerializer: a functor.  operator()(OUTPUT*, const value_type&)
    template <typename ValueSerializer, typename OUTPUT>
    bool serialize(ValueSerializer serializer, OUTPUT *fp)
    {
        if (!write_metadata(fp))
            return false;
        for (const_ne_iterator it = ne_begin(); it != ne_end(); ++it)
            if (!serializer(fp, *it))
                return false;
        return true;
    }

    // ValueSerializer: a functor.  operator()(INPUT*, value_type*)
    template <typename ValueSerializer, typename INPUT>
    bool unserialize(ValueSerializer serializer, INPUT *fp)
    {
        clear();
        if (!read_metadata(fp))
            return false;
        for (ne_iterator it = ne_begin(); it != ne_end(); ++it)
            if (!serializer(fp, &*it))
                return false;
        return true;
    }

    // Comparisons.  Note the comparisons are pretty arbitrary: we
    // compare values of the first index that isn't equal (using default
    // value for empty buckets).
    bool operator==(const sparsetable& x) const
    {
        return (_table_size == x._table_size &&
                _num_buckets == x._num_buckets &&
                _first_group == x._first_group);
    }

    bool operator<(const sparsetable& x) const
    {
        return std::lexicographical_compare(begin(), end(), x.begin(), x.end());
    }
    bool operator!=(const sparsetable& x) const { return !(*this == x); }
    bool operator<=(const sparsetable& x) const { return !(x < *this); }
    bool operator>(const sparsetable& x)  const { return x < *this; }
    bool operator>=(const sparsetable& x) const { return !(*this < x); }


private:
    // The actual data
    // ---------------
    group_type *     _first_group;
    group_type *     _last_group;
    size_type        _table_size;          // how many buckets they want
    size_type        _num_buckets;         // number of non-empty buckets
    group_alloc_type _group_alloc;
    allocator_type   _alloc;
};

//  ----------------------------------------------------------------------
//                  S P A R S E _ H A S H T A B L E
//  ----------------------------------------------------------------------
// Hashtable class, used to implement the hashed associative containers
// hash_set and hash_map.
//
// Value: what is stored in the table (each bucket is a Value).
// Key: something in a 1-to-1 correspondence to a Value, that can be used
//      to search for a Value in the table (find() takes a Key).
// HashFcn: Takes a Key and returns an integer, the more unique the better.
// ExtractKey: given a Value, returns the unique Key associated with it.
//             Must inherit from unary_function, or at least have a
//             result_type enum indicating the return type of operator().
// EqualKey: Given two Keys, says whether they are the same (that is,
//           if they are both associated with the same Value).
// Alloc: STL allocator to use to allocate memory.
//
//  ----------------------------------------------------------------------

// The probing method
// ------------------
// Linear probing
// #define JUMP_(key, num_probes)    ( 1 )
// Quadratic probing
#define JUMP_(key, num_probes)    ( num_probes )


// -------------------------------------------------------------------
// -------------------------------------------------------------------
template <class Value, class Key, class HashFcn,
          class ExtractKey, class SetKey, class EqualKey, class Alloc>
class sparse_hashtable
{
public:
    typedef Key                                        key_type;
    typedef Value                                      value_type;
    typedef HashFcn                                    hasher; // user provided or spp_hash<Key>
    typedef EqualKey                                   key_equal;
    typedef Alloc                                      allocator_type;

    typedef typename allocator_type::size_type         size_type;
    typedef typename allocator_type::difference_type   difference_type;
    typedef value_type&                                reference;
    typedef const value_type&                          const_reference;
    typedef value_type*                                pointer;
    typedef const value_type*                          const_pointer;

    // Table is the main storage class.
    typedef sparsetable<value_type, allocator_type>   Table;
    typedef typename Table::ne_iterator               ne_it;
    typedef typename Table::const_ne_iterator         cne_it;
    typedef typename Table::destructive_iterator      dest_it;
    typedef typename Table::ColIterator               ColIterator;

    typedef ne_it                                     iterator;
    typedef cne_it                                    const_iterator;
    typedef dest_it                                   destructive_iterator;

    // These come from tr1.  For us they're the same as regular iterators.
    // -------------------------------------------------------------------
    typedef iterator                                  local_iterator;
    typedef const_iterator                            const_local_iterator;

    // How full we let the table get before we resize
    // ----------------------------------------------
    static const int HT_OCCUPANCY_PCT; // = 80 (out of 100);

    // How empty we let the table get before we resize lower, by default.
    // (0.0 means never resize lower.)
    // It should be less than OCCUPANCY_PCT / 2 or we thrash resizing
    // ------------------------------------------------------------------
    static const int HT_EMPTY_PCT; // = 0.4 * HT_OCCUPANCY_PCT;

    // Minimum size we're willing to let hashtables be.
    // Must be a power of two, and at least 4.
    // Note, however, that for a given hashtable, the initial size is a
    // function of the first constructor arg, and may be >HT_MIN_BUCKETS.
    // ------------------------------------------------------------------
    static const size_type HT_MIN_BUCKETS = 4;

    // By default, if you don't specify a hashtable size at
    // construction-time, we use this size.  Must be a power of two, and
    // at least HT_MIN_BUCKETS.
    // -----------------------------------------------------------------
    static const size_type HT_DEFAULT_STARTING_BUCKETS = 32;

    // iterators
    // ---------
    iterator       begin()        { return _mk_iterator(table.ne_begin());  }
    iterator       end()          { return _mk_iterator(table.ne_end());    }
    const_iterator begin() const  { return _mk_const_iterator(table.ne_cbegin()); }
    const_iterator end() const    { return _mk_const_iterator(table.ne_cend());   }
    const_iterator cbegin() const { return _mk_const_iterator(table.ne_cbegin()); }
    const_iterator cend() const   { return _mk_const_iterator(table.ne_cend());   }

    // These come from tr1 unordered_map.  They iterate over 'bucket' n.
    // For sparsehashtable, we could consider each 'group' to be a bucket,
    // I guess, but I don't really see the point.  We'll just consider
    // bucket n to be the n-th element of the sparsetable, if it's occupied,
    // or some empty element, otherwise.
    // ---------------------------------------------------------------------
    local_iterator begin(size_type i)
    {
        return _mk_iterator(table.test(i) ? table.get_iter(i) : table.ne_end());
    }

    local_iterator end(size_type i)
    {
        local_iterator it = begin(i);
        if (table.test(i))
            ++it;
        return _mk_iterator(it);
    }

    const_local_iterator begin(size_type i) const
    {
        return _mk_const_iterator(table.test(i) ? table.get_iter(i) : table.ne_cend());
    }

    const_local_iterator end(size_type i) const
    {
        const_local_iterator it = begin(i);
        if (table.test(i))
            ++it;
        return _mk_const_iterator(it);
    }

    const_local_iterator cbegin(size_type i) const { return begin(i); }
    const_local_iterator cend(size_type i)   const { return end(i); }

    // This is used when resizing
    // --------------------------
    destructive_iterator destructive_begin()       { return _mk_destructive_iterator(table.destructive_begin()); }
    destructive_iterator destructive_end()         { return _mk_destructive_iterator(table.destructive_end());   }


    // accessor functions for the things we templatize on, basically
    // -------------------------------------------------------------
    hasher hash_funct() const               { return settings; }
    key_equal key_eq() const                { return key_info; }
    allocator_type get_allocator() const    { return table.get_allocator(); }

    // Accessor function for statistics gathering.
    unsigned int num_table_copies() const { return settings.num_ht_copies(); }

private:
    // This is used as a tag for the copy constructor, saying to destroy its
    // arg We have two ways of destructively copying: with potentially growing
    // the hashtable as we copy, and without.  To make sure the outside world
    // can't do a destructive copy, we make the typename private.
    // -----------------------------------------------------------------------
    enum MoveDontCopyT {MoveDontCopy, MoveDontGrow};

    // creating iterators from sparsetable::ne_iterators
    // -------------------------------------------------
    iterator             _mk_iterator(ne_it it) const               { return it; }
    const_iterator       _mk_const_iterator(cne_it it) const        { return it; }
    destructive_iterator _mk_destructive_iterator(dest_it it) const { return it; }

public:
    size_type size() const              { return table.num_nonempty(); }
    size_type max_size() const          { return table.max_size(); }
    bool empty() const                  { return size() == 0; }
    size_type bucket_count() const      { return table.size(); }
    size_type max_bucket_count() const  { return max_size(); }
    // These are tr1 methods.  Their idea of 'bucket' doesn't map well to
    // what we do.  We just say every bucket has 0 or 1 items in it.
    size_type bucket_size(size_type i) const
    {
        return (size_type)(begin(i) == end(i) ? 0 : 1);
    }

private:
    // Because of the above, size_type(-1) is never legal; use it for errors
    // ---------------------------------------------------------------------
    static const size_type ILLEGAL_BUCKET = size_type(-1);

    // Used after a string of deletes.  Returns true if we actually shrunk.
    // TODO(csilvers): take a delta so we can take into account inserts
    // done after shrinking.  Maybe make part of the Settings class?
    // --------------------------------------------------------------------
    bool _maybe_shrink()
    {
        assert((bucket_count() & (bucket_count()-1)) == 0); // is a power of two
        assert(bucket_count() >= HT_MIN_BUCKETS);
        bool retval = false;

        // If you construct a hashtable with < HT_DEFAULT_STARTING_BUCKETS,
        // we'll never shrink until you get relatively big, and we'll never
        // shrink below HT_DEFAULT_STARTING_BUCKETS.  Otherwise, something
        // like "dense_hash_set<int> x; x.insert(4); x.erase(4);" will
        // shrink us down to HT_MIN_BUCKETS buckets, which is too small.
        // ---------------------------------------------------------------
        const size_type num_remain = table.num_nonempty();
        const size_type shrink_threshold = settings.shrink_threshold();
        if (shrink_threshold > 0 && num_remain < shrink_threshold &&
            bucket_count() > HT_DEFAULT_STARTING_BUCKETS)
        {
            const float shrink_factor = settings.shrink_factor();
            size_type sz = (size_type)(bucket_count() / 2);    // find how much we should shrink
            while (sz > HT_DEFAULT_STARTING_BUCKETS &&
                   num_remain < static_cast<size_type>(sz * shrink_factor))
            {
                sz /= 2;                            // stay a power of 2
            }
            sparse_hashtable tmp(MoveDontCopy, *this, sz);
            swap(tmp);                            // now we are tmp
            retval = true;
        }
        settings.set_consider_shrink(false);   // because we just considered it
        return retval;
    }

    // We'll let you resize a hashtable -- though this makes us copy all!
    // When you resize, you say, "make it big enough for this many more elements"
    // Returns true if we actually resized, false if size was already ok.
    // --------------------------------------------------------------------------
    bool _resize_delta(size_type delta)
    {
        bool did_resize = false;
        if (settings.consider_shrink())
        {
            // see if lots of deletes happened
            if (_maybe_shrink())
                did_resize = true;
        }
        if (table.num_nonempty() >=
            (std::numeric_limits<size_type>::max)() - delta)
        {
            throw_exception(std::length_error("resize overflow"));
        }

        size_type num_occupied = (size_type)(table.num_nonempty() + num_deleted);

        if (bucket_count() >= HT_MIN_BUCKETS &&
             (num_occupied + delta) <= settings.enlarge_threshold())
            return did_resize;                       // we're ok as we are

        // Sometimes, we need to resize just to get rid of all the
        // "deleted" buckets that are clogging up the hashtable.  So when
        // deciding whether to resize, count the deleted buckets (which
        // are currently taking up room).
        // -------------------------------------------------------------
        const size_type needed_size =
                  settings.min_buckets((size_type)(num_occupied + delta), (size_type)0);

        if (needed_size <= bucket_count())      // we have enough buckets
            return did_resize;

        size_type resize_to = settings.min_buckets((size_type)(num_occupied + delta), bucket_count());

        if (resize_to < needed_size &&    // may double resize_to
            resize_to < (std::numeric_limits<size_type>::max)() / 2)
        {
            // This situation means that we have enough deleted elements,
            // that once we purge them, we won't actually have needed to
            // grow.  But we may want to grow anyway: if we just purge one
            // element, say, we'll have to grow anyway next time we
            // insert.  Might as well grow now, since we're already going
            // through the trouble of copying (in order to purge the
            // deleted elements).
            const size_type target =
                static_cast<size_type>(settings.shrink_size((size_type)(resize_to*2)));
            if (table.num_nonempty() + delta >= target)
            {
                // Good, we won't be below the shrink threshhold even if we double.
                resize_to *= 2;
            }
        }

        sparse_hashtable tmp(MoveDontCopy, *this, resize_to);
        swap(tmp);                             // now we are tmp
        return true;
    }

    // Used to actually do the rehashing when we grow/shrink a hashtable
    // -----------------------------------------------------------------
    void _copy_from(const sparse_hashtable &ht, size_type min_buckets_wanted)
    {
        clear();            // clear table, set num_deleted to 0

        // If we need to change the size of our table, do it now
        const size_type resize_to = settings.min_buckets(ht.size(), min_buckets_wanted);

        if (resize_to > bucket_count())
        {
            // we don't have enough buckets
            table.resize(resize_to);               // sets the number of buckets
            settings.reset_thresholds(bucket_count());
        }

        // We use a normal iterator to get bcks from ht
        // We could use insert() here, but since we know there are
        // no duplicates, we can be more efficient
        assert((bucket_count() & (bucket_count()-1)) == 0);      // a power of two
        for (const_iterator it = ht.begin(); it != ht.end(); ++it)
        {
            size_type num_probes = 0;              // how many times we've probed
            size_type bucknum;
            const size_type bucket_count_minus_one = bucket_count() - 1;
            for (bucknum = hash(get_key(*it)) & bucket_count_minus_one;
                 table.test(bucknum);                                   // table.test() OK since no erase()
                 bucknum = (bucknum + JUMP_(key, num_probes)) & bucket_count_minus_one)
            {
                ++num_probes;
                assert(num_probes < bucket_count()
                       && "Hashtable is full: an error in key_equal<> or hash<>");
            }
            table.set(bucknum, *it);               // copies the value to here
        }
        settings.inc_num_ht_copies();
    }

    // Implementation is like _copy_from, but it destroys the table of the
    // "from" guy by freeing sparsetable memory as we iterate.  This is
    // useful in resizing, since we're throwing away the "from" guy anyway.
    // --------------------------------------------------------------------
    void _move_from(MoveDontCopyT mover, sparse_hashtable &ht,
                   size_type min_buckets_wanted)
    {
        clear();

        // If we need to change the size of our table, do it now
        size_type resize_to;
        if (mover == MoveDontGrow)
            resize_to = ht.bucket_count();       // keep same size as old ht
        else                                     // MoveDontCopy
            resize_to = settings.min_buckets(ht.size(), min_buckets_wanted);
        if (resize_to > bucket_count())
        {
            // we don't have enough buckets
            table.resize(resize_to);               // sets the number of buckets
            settings.reset_thresholds(bucket_count());
        }

        // We use a normal iterator to get bcks from ht
        // We could use insert() here, but since we know there are
        // no duplicates, we can be more efficient
        assert((bucket_count() & (bucket_count()-1)) == 0);      // a power of two
        const size_type bucket_count_minus_one = (const size_type)(bucket_count() - 1);

        // THIS IS THE MAJOR LINE THAT DIFFERS FROM COPY_FROM():
        for (destructive_iterator it = ht.destructive_begin();
              it != ht.destructive_end(); ++it)
        {
            size_type num_probes = 0;
            size_type bucknum;
            for (bucknum = hash(get_key(*it)) & bucket_count_minus_one;
                 table.test(bucknum);                          // table.test() OK since no erase()
                 bucknum = (size_type)((bucknum + JUMP_(key, num_probes)) & (bucket_count()-1)))
            {
                ++num_probes;
                assert(num_probes < bucket_count()
                       && "Hashtable is full: an error in key_equal<> or hash<>");
            }
            table.move(bucknum, *it);    // moves the value to here
        }
        settings.inc_num_ht_copies();
    }


    // Required by the spec for hashed associative container
public:
    // Though the docs say this should be num_buckets, I think it's much
    // more useful as num_elements.  As a special feature, calling with
    // req_elements==0 will cause us to shrink if we can, saving space.
    // -----------------------------------------------------------------
    void resize(size_type req_elements)
    {
        // resize to this or larger
        if (settings.consider_shrink() || req_elements == 0)
            _maybe_shrink();
        if (req_elements > table.num_nonempty())    // we only grow
            _resize_delta((size_type)(req_elements - table.num_nonempty()));
    }

    // Get and change the value of shrink_factor and enlarge_factor.  The
    // description at the beginning of this file explains how to choose
    // the values.  Setting the shrink parameter to 0.0 ensures that the
    // table never shrinks.
    // ------------------------------------------------------------------
    void get_resizing_parameters(float* shrink, float* grow) const
    {
        *shrink = settings.shrink_factor();
        *grow = settings.enlarge_factor();
    }

    float get_shrink_factor() const  { return settings.shrink_factor(); }
    float get_enlarge_factor() const { return settings.enlarge_factor(); }

    void set_resizing_parameters(float shrink, float grow) 
    {
        settings.set_resizing_parameters(shrink, grow);
        settings.reset_thresholds(bucket_count());
    }

    void set_shrink_factor(float shrink)
    {
        set_resizing_parameters(shrink, get_enlarge_factor());
    }

    void set_enlarge_factor(float grow)
    {
        set_resizing_parameters(get_shrink_factor(), grow);
    }

    // CONSTRUCTORS -- as required by the specs, we take a size,
    // but also let you specify a hashfunction, key comparator,
    // and key extractor.  We also define a copy constructor and =.
    // DESTRUCTOR -- the default is fine, surprisingly.
    // ------------------------------------------------------------
    explicit sparse_hashtable(size_type expected_max_items_in_table = 0,
                              const HashFcn& hf = HashFcn(),
                              const EqualKey& eql = EqualKey(),
                              const ExtractKey& ext = ExtractKey(),
                              const SetKey& set = SetKey(),
                              const allocator_type& alloc = allocator_type())
        : settings(hf),
          key_info(ext, set, eql),
          num_deleted(0),
          table((expected_max_items_in_table == 0
                 ? HT_DEFAULT_STARTING_BUCKETS
                 : settings.min_buckets(expected_max_items_in_table, 0)),
                alloc)
    {
        settings.reset_thresholds(bucket_count());
    }

    // As a convenience for resize(), we allow an optional second argument
    // which lets you make this new hashtable a different size than ht.
    // We also provide a mechanism of saying you want to "move" the ht argument
    // into us instead of copying.
    // ------------------------------------------------------------------------
    sparse_hashtable(const sparse_hashtable& ht,
                     size_type min_buckets_wanted = HT_DEFAULT_STARTING_BUCKETS)
        : settings(ht.settings),
          key_info(ht.key_info),
          num_deleted(0),
          table(0)
    {
        settings.reset_thresholds(bucket_count());
        _copy_from(ht, min_buckets_wanted);
    }

#if !defined(SPP_NO_CXX11_RVALUE_REFERENCES)

    sparse_hashtable(sparse_hashtable&& o, const allocator_type& alloc = allocator_type()) :
        settings(o.settings),
        key_info(o.key_info),
        num_deleted(0),
        table(HT_DEFAULT_STARTING_BUCKETS, alloc)
    {
        settings.reset_thresholds(bucket_count());
        this->swap(o);
    }

    sparse_hashtable& operator=(sparse_hashtable&& o)
    {
        this->swap(o);
        return *this;
    }
#endif

    sparse_hashtable(MoveDontCopyT mover,
                     sparse_hashtable& ht,
                     size_type min_buckets_wanted = HT_DEFAULT_STARTING_BUCKETS)
        : settings(ht.settings),
          key_info(ht.key_info),
          num_deleted(0),
          table(min_buckets_wanted, ht.table.get_allocator())
          //table(min_buckets_wanted)
    {
        settings.reset_thresholds(bucket_count());
        _move_from(mover, ht, min_buckets_wanted);
    }

    sparse_hashtable& operator=(const sparse_hashtable& ht)
    {
        if (&ht == this)
            return *this;        // don't copy onto ourselves
        settings = ht.settings;
        key_info = ht.key_info;
        num_deleted = ht.num_deleted;

        // _copy_from() calls clear and sets num_deleted to 0 too
        _copy_from(ht, HT_MIN_BUCKETS);

        // we purposefully don't copy the allocator, which may not be copyable
        return *this;
    }

    // Many STL algorithms use swap instead of copy constructors
    void swap(sparse_hashtable& ht)
    {
        using std::swap;

        swap(settings, ht.settings);
        swap(key_info, ht.key_info);
        swap(num_deleted, ht.num_deleted);
        table.swap(ht.table);
        settings.reset_thresholds(bucket_count());  // also resets consider_shrink
        ht.settings.reset_thresholds(ht.bucket_count());
        // we purposefully don't swap the allocator, which may not be swap-able
    }

    // It's always nice to be able to clear a table without deallocating it
    void clear()
    {
        if (!empty() || num_deleted != 0)
        {
            table.clear();
            table = Table(HT_DEFAULT_STARTING_BUCKETS, table.get_allocator());
        }
        settings.reset_thresholds(bucket_count());
        num_deleted = 0;
    }

    // LOOKUP ROUTINES
private:

    enum pos_type { pt_empty = 0, pt_erased, pt_full };
    // -------------------------------------------------------------------
    class Position
    {
    public:

        Position() : _t(pt_empty) {}
        Position(pos_type t, size_type idx) : _t(t), _idx(idx) {}

        pos_type  _t;
        size_type _idx;
    };

    // Returns a pair:
    //   - 'first' is a code, 2 if key already present, 0 or 1 otherwise.
    //   - 'second' is a position, where the key should go
    // Note: because of deletions where-to-insert is not trivial: it's the
    // first deleted bucket we see, as long as we don't find the key later
    // -------------------------------------------------------------------
    Position _find_position(const key_type &key) const
    {
        size_type num_probes = 0;                    // how many times we've probed
        const size_type bucket_count_minus_one = (const size_type)(bucket_count() - 1);
        size_type bucknum = hash(key) & bucket_count_minus_one;
        Position pos;

        while (1)
        {
            // probe until something happens
            // -----------------------------
            typename Table::GrpPos grp_pos(table, bucknum);

            if (!grp_pos.test_strict())
            {
                // bucket is empty => key not present
                return pos._t ? pos : Position(pt_empty, bucknum);
            }
            else if (grp_pos.test())
            {
                reference ref(grp_pos.unsafe_get());

                if (equals(key, get_key(ref)))
                    return Position(pt_full, bucknum);
            }
            else if (pos._t == pt_empty)
            {
                // first erased position
                pos._t   = pt_erased;
                pos._idx = bucknum;
            }

            ++num_probes;                        // we're doing another probe
            bucknum = (size_type)((bucknum + JUMP_(key, num_probes)) & bucket_count_minus_one);
            assert(num_probes < bucket_count()
                   && "Hashtable is full: an error in key_equal<> or hash<>");
        }
    }

public:
    // I hate to duplicate find() like that, but it is
    // significantly faster to not have the intermediate pair
    // ------------------------------------------------------------------
    iterator find(const key_type& key)
    {
        size_type num_probes = 0;              // how many times we've probed
        const size_type bucket_count_minus_one = bucket_count() - 1;
        size_type bucknum = hash(key) & bucket_count_minus_one;

        while (1)                        // probe until something happens
        {
            typename Table::GrpPos grp_pos(table, bucknum);

            if (!grp_pos.test_strict())
                return end();            // bucket is empty
            if (grp_pos.test())
            {
                reference ref(grp_pos.unsafe_get());

                if (equals(key, get_key(ref)))
                    return grp_pos.get_iter(ref);
            }
            ++num_probes;                        // we're doing another probe
            bucknum = (bucknum + JUMP_(key, num_probes)) & bucket_count_minus_one;
            assert(num_probes < bucket_count()
                   && "Hashtable is full: an error in key_equal<> or hash<>");
        }
    }

    // Wish I could avoid the duplicate find() const and non-const.
    // ------------------------------------------------------------
    const_iterator find(const key_type& key) const
    {
        size_type num_probes = 0;              // how many times we've probed
        const size_type bucket_count_minus_one = bucket_count() - 1;
        size_type bucknum = hash(key) & bucket_count_minus_one;

        while (1)                        // probe until something happens
        {
            typename Table::GrpPos grp_pos(table, bucknum);

            if (!grp_pos.test_strict())
                return end();            // bucket is empty
            else if (grp_pos.test())
            {
                reference ref(grp_pos.unsafe_get());

                if (equals(key, get_key(ref)))
                    return _mk_const_iterator(table.get_iter(bucknum, &ref));
            }
            ++num_probes;                        // we're doing another probe
            bucknum = (bucknum + JUMP_(key, num_probes)) & bucket_count_minus_one;
            assert(num_probes < bucket_count()
                   && "Hashtable is full: an error in key_equal<> or hash<>");
        }
    }

    // This is a tr1 method: the bucket a given key is in, or what bucket
    // it would be put in, if it were to be inserted.  Shrug.
    // ------------------------------------------------------------------
    size_type bucket(const key_type& key) const
    {
        Position pos = _find_position(key);
        return pos._idx;
    }

    // Counts how many elements have key key.  For maps, it's either 0 or 1.
    // ---------------------------------------------------------------------
    size_type count(const key_type &key) const
    {
        Position pos = _find_position(key);
        return (size_type)(pos._t == pt_full ? 1 : 0);
    }

    // Likewise, equal_range doesn't really make sense for us.  Oh well.
    // -----------------------------------------------------------------
    std::pair<iterator,iterator> equal_range(const key_type& key)
    {
        iterator pos = find(key);      // either an iterator or end
        if (pos == end())
            return std::pair<iterator,iterator>(pos, pos);
        else
        {
            const iterator startpos = pos++;
            return std::pair<iterator,iterator>(startpos, pos);
        }
    }

    std::pair<const_iterator,const_iterator> equal_range(const key_type& key) const
    {
        const_iterator pos = find(key);      // either an iterator or end
        if (pos == end())
            return std::pair<const_iterator,const_iterator>(pos, pos);
        else
        {
            const const_iterator startpos = pos++;
            return std::pair<const_iterator,const_iterator>(startpos, pos);
        }
    }


    // INSERTION ROUTINES
private:
    // Private method used by insert_noresize and find_or_insert.
    template <class T>
    reference _insert_at(T& obj, size_type pos, bool erased)
    {
        if (size() >= max_size())
        {
            throw_exception(std::length_error("insert overflow"));
        }
        if (erased)
        {
            assert(num_deleted);
            --num_deleted;
        }
        return table.set(pos, obj);
    }

    // If you know *this is big enough to hold obj, use this routine
    template <class T>
    std::pair<iterator, bool> _insert_noresize(T& obj)
    {
        Position pos = _find_position(get_key(obj));
        bool already_there = (pos._t == pt_full);

        if (!already_there)
        {
            reference ref(_insert_at(obj, pos._idx, pos._t == pt_erased));
            return std::pair<iterator, bool>(_mk_iterator(table.get_iter(pos._idx, &ref)), true);
        }
        return std::pair<iterator,bool>(_mk_iterator(table.get_iter(pos._idx)), false);
    }

    // Specializations of insert(it, it) depending on the power of the iterator:
    // (1) Iterator supports operator-, resize before inserting
    template <class ForwardIterator>
    void _insert(ForwardIterator f, ForwardIterator l, std::forward_iterator_tag /*unused*/)
    {
        int64_t dist = std::distance(f, l);
        if (dist < 0 ||  static_cast<size_t>(dist) >= (std::numeric_limits<size_type>::max)())
            throw_exception(std::length_error("insert-range overflow"));

        _resize_delta(static_cast<size_type>(dist));

        for (; dist > 0; --dist, ++f)
            _insert_noresize(*f);
    }

    // (2) Arbitrary iterator, can't tell how much to resize
    template <class InputIterator>
    void _insert(InputIterator f, InputIterator l, std::input_iterator_tag /*unused*/)
    {
        for (; f != l; ++f)
            _insert(*f);
    }

public:

#if !defined(SPP_NO_CXX11_VARIADIC_TEMPLATES)
    template <class... Args>
    std::pair<iterator, bool> emplace(Args&&... args)
    {
        _resize_delta(1);
        value_type obj(std::forward<Args>(args)...);
        return _insert_noresize(obj);
    }
#endif

    // This is the normal insert routine, used by the outside world
    std::pair<iterator, bool> insert(const_reference obj)
    {
        _resize_delta(1);                      // adding an object, grow if need be
        return _insert_noresize(obj);
    }

#if !defined(SPP_NO_CXX11_RVALUE_REFERENCES)
    template< class P >
    std::pair<iterator, bool> insert(P &&obj)
    {
        _resize_delta(1);                      // adding an object, grow if need be
        value_type val(std::forward<P>(obj));
        return _insert_noresize(val);
    }
#endif

    // When inserting a lot at a time, we specialize on the type of iterator
    template <class InputIterator>
    void insert(InputIterator f, InputIterator l)
    {
        // specializes on iterator type
        _insert(f, l,
               typename std::iterator_traits<InputIterator>::iterator_category());
    }

    // DefaultValue is a functor that takes a key and returns a value_type
    // representing the default value to be inserted if none is found.
#if !defined(SPP_NO_CXX11_VARIADIC_TEMPLATES)
    template <class DefaultValue, class KT>
    value_type& find_or_insert(KT&& key)
#else
    template <class DefaultValue>
    value_type& find_or_insert(const key_type& key)
#endif
    {
        size_type num_probes = 0;              // how many times we've probed
        const size_type bucket_count_minus_one = bucket_count() - 1;
        size_type bucknum = hash(key) & bucket_count_minus_one;
        DefaultValue default_value;
        size_type erased_pos = 0;
        bool erased = false;

        while (1)                        // probe until something happens
        {
            typename Table::GrpPos grp_pos(table, bucknum);

            if (!grp_pos.test_strict())
            {
                // not found
#if !defined(SPP_NO_CXX11_VARIADIC_TEMPLATES)
                auto&& def(default_value(std::forward<KT>(key)));
#else
                value_type def(default_value(key));
#endif                
                if (_resize_delta(1))
                {
                    // needed to rehash to make room
                    // Since we resized, we can't use pos, so recalculate where to insert.
                    return *(_insert_noresize(def).first);
                }
                else
                {
                    // no need to rehash, insert right here
                    return _insert_at(def, erased ? erased_pos : bucknum, erased);
                }
            }
            if (grp_pos.test())
            {
                reference ref(grp_pos.unsafe_get());

                if (equals(key, get_key(ref)))
                    return ref;
            }
            else if (!erased)
            {
                // first erased position
                erased_pos = bucknum;
                erased = true;
            }

            ++num_probes;                        // we're doing another probe
            bucknum = (bucknum + JUMP_(key, num_probes)) & bucket_count_minus_one;
            assert(num_probes < bucket_count()
                   && "Hashtable is full: an error in key_equal<> or hash<>");
        }
    }

    size_type erase(const key_type& key)
    {
        size_type num_probes = 0;              // how many times we've probed
        const size_type bucket_count_minus_one = bucket_count() - 1;
        size_type bucknum = hash(key) & bucket_count_minus_one;

        while (1)                        // probe until something happens
        {
            typename Table::GrpPos grp_pos(table, bucknum);

            if (!grp_pos.test_strict())
                return 0;            // bucket is empty, we deleted nothing
            if (grp_pos.test())
            {
                reference ref(grp_pos.unsafe_get());

                if (equals(key, get_key(ref)))
                {
                    grp_pos.erase(table);
                    ++num_deleted;
                    settings.set_consider_shrink(true); // will think about shrink after next insert
                    return 1;                           // because we deleted one thing
                }
            }
            ++num_probes;                        // we're doing another probe
            bucknum = (bucknum + JUMP_(key, num_probes)) & bucket_count_minus_one;
            assert(num_probes < bucket_count()
                   && "Hashtable is full: an error in key_equal<> or hash<>");
        }
    }

    const_iterator erase(const_iterator pos)
    {
        if (pos == cend())
            return cend();                 // sanity check

        const_iterator nextpos = table.erase(pos);
        ++num_deleted;
        settings.set_consider_shrink(true);
        return nextpos;
    }

    const_iterator erase(const_iterator f, const_iterator l)
    {
        if (f == cend())
            return cend();                // sanity check

        size_type num_before = table.num_nonempty();
        const_iterator nextpos = table.erase(f, l);
        num_deleted += num_before - table.num_nonempty();
        settings.set_consider_shrink(true);
        return nextpos;
    }

    // Deleted key routines - just to keep google test framework happy
    // we don't actually use the deleted key
    // ---------------------------------------------------------------
    void set_deleted_key(const key_type&)
    {
    }

    void clear_deleted_key()
    {
    }

    bool operator==(const sparse_hashtable& ht) const
    {
        if (this == &ht)
            return true;

        if (size() != ht.size())
            return false;

        for (const_iterator it = begin(); it != end(); ++it)
        {
            const_iterator it2 = ht.find(get_key(*it));
            if ((it2 == ht.end()) || (*it != *it2))
                return false;
        }

        return true;
    }

    bool operator!=(const sparse_hashtable& ht) const
    {
        return !(*this == ht);
    }


    // I/O
    // We support reading and writing hashtables to disk.  NOTE that
    // this only stores the hashtable metadata, not the stuff you've
    // actually put in the hashtable!  Alas, since I don't know how to
    // write a hasher or key_equal, you have to make sure everything
    // but the table is the same.  We compact before writing.
    //
    // The OUTPUT type needs to support a Write() operation. File and
    // OutputBuffer are appropriate types to pass in.
    //
    // The INPUT type needs to support a Read() operation. File and
    // InputBuffer are appropriate types to pass in.
    // -------------------------------------------------------------
    template <typename OUTPUT>
    bool write_metadata(OUTPUT *fp)
    {
        return table.write_metadata(fp);
    }

    template <typename INPUT>
    bool read_metadata(INPUT *fp)
    {
        num_deleted = 0;            // since we got rid before writing
        const bool result = table.read_metadata(fp);
        settings.reset_thresholds(bucket_count());
        return result;
    }

    // Only meaningful if value_type is a POD.
    template <typename OUTPUT>
    bool write_nopointer_data(OUTPUT *fp)
    {
        return table.write_nopointer_data(fp);
    }

    // Only meaningful if value_type is a POD.
    template <typename INPUT>
    bool read_nopointer_data(INPUT *fp)
    {
        return table.read_nopointer_data(fp);
    }

    // INPUT and OUTPUT must be either a FILE, *or* a C++ stream
    //    (istream, ostream, etc) *or* a class providing
    //    Read(void*, size_t) and Write(const void*, size_t)
    //    (respectively), which writes a buffer into a stream
    //    (which the INPUT/OUTPUT instance presumably owns).

    typedef sparsehash_internal::pod_serializer<value_type> NopointerSerializer;

    // ValueSerializer: a functor.  operator()(OUTPUT*, const value_type&)
    template <typename ValueSerializer, typename OUTPUT>
    bool serialize(ValueSerializer serializer, OUTPUT *fp)
    {
        return table.serialize(serializer, fp);
    }

    // ValueSerializer: a functor.  operator()(INPUT*, value_type*)
    template <typename ValueSerializer, typename INPUT>
    bool unserialize(ValueSerializer serializer, INPUT *fp)
    {
        num_deleted = 0;            // since we got rid before writing
        const bool result = table.unserialize(serializer, fp);
        settings.reset_thresholds(bucket_count());
        return result;
    }

private:

    // Package templated functors with the other types to eliminate memory
    // needed for storing these zero-size operators.  Since ExtractKey and
    // hasher's operator() might have the same function signature, they
    // must be packaged in different classes.
    // -------------------------------------------------------------------------
    struct Settings :
        sparsehash_internal::sh_hashtable_settings<key_type, hasher,
                                                   size_type, HT_MIN_BUCKETS>
    {
        explicit Settings(const hasher& hf)
            : sparsehash_internal::sh_hashtable_settings<key_type, hasher, size_type,
              HT_MIN_BUCKETS>
              (hf, HT_OCCUPANCY_PCT / 100.0f, HT_EMPTY_PCT / 100.0f) {}
    };

    // KeyInfo stores delete key and packages zero-size functors:
    // ExtractKey and SetKey.
     // ---------------------------------------------------------
    class KeyInfo : public ExtractKey, public SetKey, public EqualKey
    {
    public:
        KeyInfo(const ExtractKey& ek, const SetKey& sk, const EqualKey& eq)
            : ExtractKey(ek), SetKey(sk), EqualKey(eq)
        {
        }

        // We want to return the exact same type as ExtractKey: Key or const Key&
        typename ExtractKey::result_type get_key(const_reference v) const
        {
            return ExtractKey::operator()(v);
        }

        bool equals(const key_type& a, const key_type& b) const
        {
            return EqualKey::operator()(a, b);
        }
    };

    // Utility functions to access the templated operators
    size_t hash(const key_type& v) const
    {
        return settings.hash(v);
    }

    bool equals(const key_type& a, const key_type& b) const
    {
        return key_info.equals(a, b);
    }

    typename ExtractKey::result_type get_key(const_reference v) const
    {
        return key_info.get_key(v);
    }

private:
    // Actual data
    // -----------
    Settings  settings;
    KeyInfo   key_info;
    size_type num_deleted;
    Table     table;         // holds num_buckets and num_elements too
};

#undef JUMP_

// -----------------------------------------------------------------------------
template <class V, class K, class HF, class ExK, class SetK, class EqK, class A>
const typename sparse_hashtable<V,K,HF,ExK,SetK,EqK,A>::size_type
sparse_hashtable<V,K,HF,ExK,SetK,EqK,A>::ILLEGAL_BUCKET;

// How full we let the table get before we resize.  Knuth says .8 is
// good -- higher causes us to probe too much, though saves memory
// -----------------------------------------------------------------------------
template <class V, class K, class HF, class ExK, class SetK, class EqK, class A>
const int sparse_hashtable<V,K,HF,ExK,SetK,EqK,A>::HT_OCCUPANCY_PCT = 50;

// How empty we let the table get before we resize lower.
// It should be less than OCCUPANCY_PCT / 2 or we thrash resizing
// -----------------------------------------------------------------------------
template <class V, class K, class HF, class ExK, class SetK, class EqK, class A>
const int sparse_hashtable<V,K,HF,ExK,SetK,EqK,A>::HT_EMPTY_PCT
= static_cast<int>(0.4 *
                   sparse_hashtable<V,K,HF,ExK,SetK,EqK,A>::HT_OCCUPANCY_PCT);


//  ----------------------------------------------------------------------
//                   S P A R S E _ H A S H _ M A P
//  ----------------------------------------------------------------------
template <class Key, class T,
          class HashFcn  = spp_hash<Key>,
          class EqualKey = std::equal_to<Key>,
          class Alloc    = SPP_DEFAULT_ALLOCATOR<std::pair<const Key, T> > >
class sparse_hash_map
{
public:
    typedef typename std::pair<const Key, T> value_type;

private:
    // Apparently select1st is not stl-standard, so we define our own
    struct SelectKey
    {
        typedef const Key& result_type;

        inline const Key& operator()(const value_type& p) const
        {
            return p.first;
        }
    };

    struct SetKey
    {
        inline void operator()(value_type* value, const Key& new_key) const
        {
            *const_cast<Key*>(&value->first) = new_key;
        }
    };

    // For operator[].
    struct DefaultValue
    {
#if !defined(SPP_NO_CXX11_VARIADIC_TEMPLATES)
        template <class KT>
        inline value_type operator()(KT&& key)  const
        {
            return { std::forward<KT>(key), T() };
        }
#else
        inline value_type operator()(const Key& key)  const
        {
            return std::make_pair(key, T());
        }
#endif
    };

    // The actual data
    typedef sparse_hashtable<value_type, Key, HashFcn, SelectKey,
                             SetKey, EqualKey, Alloc> ht;

public:
    typedef typename ht::key_type             key_type;
    typedef T                                 data_type;
    typedef T                                 mapped_type;
    typedef typename ht::hasher               hasher;
    typedef typename ht::key_equal            key_equal;
    typedef Alloc                             allocator_type;

    typedef typename ht::size_type            size_type;
    typedef typename ht::difference_type      difference_type;
    typedef typename ht::pointer              pointer;
    typedef typename ht::const_pointer        const_pointer;
    typedef typename ht::reference            reference;
    typedef typename ht::const_reference      const_reference;

    typedef typename ht::iterator             iterator;
    typedef typename ht::const_iterator       const_iterator;
    typedef typename ht::local_iterator       local_iterator;
    typedef typename ht::const_local_iterator const_local_iterator;

    // Iterator functions
    iterator       begin()                         { return rep.begin(); }
    iterator       end()                           { return rep.end(); }
    const_iterator begin() const                   { return rep.cbegin(); }
    const_iterator end() const                     { return rep.cend(); }
    const_iterator cbegin() const                  { return rep.cbegin(); }
    const_iterator cend() const                    { return rep.cend(); }

    // These come from tr1's unordered_map. For us, a bucket has 0 or 1 elements.
    local_iterator begin(size_type i)              { return rep.begin(i); }
    local_iterator end(size_type i)                { return rep.end(i); }
    const_local_iterator begin(size_type i) const  { return rep.begin(i); }
    const_local_iterator end(size_type i) const    { return rep.end(i); }
    const_local_iterator cbegin(size_type i) const { return rep.cbegin(i); }
    const_local_iterator cend(size_type i) const   { return rep.cend(i); }

    // Accessor functions
    // ------------------
    allocator_type get_allocator() const           { return rep.get_allocator(); }
    hasher hash_funct() const                      { return rep.hash_funct(); }
    hasher hash_function() const                   { return hash_funct(); }
    key_equal key_eq() const                       { return rep.key_eq(); }


    // Constructors
    // ------------
    explicit sparse_hash_map(size_type n = 0,
                             const hasher& hf = hasher(),
                             const key_equal& eql = key_equal(),
                             const allocator_type& alloc = allocator_type())
        : rep(n, hf, eql, SelectKey(), SetKey(), alloc)
    {
    }

    explicit sparse_hash_map(const allocator_type& alloc) :
        rep(0, hasher(), key_equal(), SelectKey(), SetKey(), alloc)
    {
    }

    sparse_hash_map(size_type n, const allocator_type& alloc) :
        rep(n, hasher(), key_equal(), SelectKey(), SetKey(), alloc)
    {
    }

    sparse_hash_map(size_type n, const hasher& hf, const allocator_type& alloc) :
        rep(n, hf, key_equal(), SelectKey(), SetKey(), alloc)
    {
    }

    template <class InputIterator>
    sparse_hash_map(InputIterator f, InputIterator l,
                    size_type n = 0,
                    const hasher& hf = hasher(),
                    const key_equal& eql = key_equal(),
                    const allocator_type& alloc = allocator_type())
        : rep(n, hf, eql, SelectKey(), SetKey(), alloc)
    {
        rep.insert(f, l);
    }

    template <class InputIterator>
    sparse_hash_map(InputIterator f, InputIterator l,
                    size_type n, const allocator_type& alloc)
        : rep(n, hasher(), key_equal(), SelectKey(), SetKey(), alloc)
    {
        rep.insert(f, l);
    }

    template <class InputIterator>
    sparse_hash_map(InputIterator f, InputIterator l,
                    size_type n, const hasher& hf, const allocator_type& alloc)
        : rep(n, hf, key_equal(), SelectKey(), SetKey(), alloc)
    {
        rep.insert(f, l);
    }

    sparse_hash_map(const sparse_hash_map &o) :
        rep(o.rep)
    {}

    sparse_hash_map(const sparse_hash_map &o,
                    const allocator_type& alloc) :
        rep(o.rep, alloc)
    {}

#if !defined(SPP_NO_CXX11_RVALUE_REFERENCES)
    sparse_hash_map(sparse_hash_map &&o) :
        rep(std::move(o.rep))
    {}

    sparse_hash_map(sparse_hash_map &&o,
                    const allocator_type& alloc) :
        rep(std::move(o.rep), alloc)
    {}

    sparse_hash_map& operator=(sparse_hash_map &&o) = default;
#endif

#if !defined(SPP_NO_CXX11_HDR_INITIALIZER_LIST)
    sparse_hash_map(std::initializer_list<value_type> init,
                    size_type n = 0,
                    const hasher& hf = hasher(),
                    const key_equal& eql = key_equal(),
                    const allocator_type& alloc = allocator_type())
        : rep(n, hf, eql, SelectKey(), SetKey(), alloc)
    {
        rep.insert(init.begin(), init.end());
    }

    sparse_hash_map(std::initializer_list<value_type> init,
                    size_type n, const allocator_type& alloc) :
        rep(n, hasher(), key_equal(), SelectKey(), SetKey(), alloc)
    {
        rep.insert(init.begin(), init.end());
    }

    sparse_hash_map(std::initializer_list<value_type> init,
                    size_type n, const hasher& hf, const allocator_type& alloc) :
        rep(n, hf, key_equal(), SelectKey(), SetKey(), alloc)
    {
        rep.insert(init.begin(), init.end());
    }

    sparse_hash_map& operator=(std::initializer_list<value_type> init)
    {
        rep.clear();
        rep.insert(init.begin(), init.end());
        return *this;
    }

    void insert(std::initializer_list<value_type> init)
    {
        rep.insert(init.begin(), init.end());
    }
#endif

    sparse_hash_map& operator=(const sparse_hash_map &o)
    {
        rep = o.rep;
        return *this;
    }

    void clear()                        { rep.clear(); }
    void swap(sparse_hash_map& hs)      { rep.swap(hs.rep); }

    // Functions concerning size
    // -------------------------
    size_type size() const              { return rep.size(); }
    size_type max_size() const          { return rep.max_size(); }
    bool empty() const                  { return rep.empty(); }
    size_type bucket_count() const      { return rep.bucket_count(); }
    size_type max_bucket_count() const  { return rep.max_bucket_count(); }

    size_type bucket_size(size_type i) const    { return rep.bucket_size(i); }
    size_type bucket(const key_type& key) const { return rep.bucket(key); }
    float     load_factor() const       { return size() * 1.0f / bucket_count(); }

    float max_load_factor() const      { return rep.get_enlarge_factor(); }
    void  max_load_factor(float grow)  { rep.set_enlarge_factor(grow); }

    float min_load_factor() const      { return rep.get_shrink_factor(); }
    void  min_load_factor(float shrink){ rep.set_shrink_factor(shrink); }

    void set_resizing_parameters(float shrink, float grow)
    {
        rep.set_resizing_parameters(shrink, grow);
    }

    void resize(size_type cnt)        { rep.resize(cnt); }
    void rehash(size_type cnt)        { resize(cnt); } // c++11 name
    void reserve(size_type cnt)       { resize(cnt); } // c++11

    // Lookup
    // ------
    iterator find(const key_type& key)                 { return rep.find(key); }
    const_iterator find(const key_type& key) const     { return rep.find(key); }
    bool contains(const key_type& key) const           { return rep.find(key) != rep.end(); }

#if !defined(SPP_NO_CXX11_VARIADIC_TEMPLATES)
    template <class KT>
    mapped_type& operator[](KT&& key)
    {
        return rep.template find_or_insert<DefaultValue>(std::forward<KT>(key)).second;
    }
#else
    mapped_type& operator[](const key_type& key)
    {
        return rep.template find_or_insert<DefaultValue>(key).second;
    }
#endif

    size_type count(const key_type& key) const         { return rep.count(key); }

    std::pair<iterator, iterator>
    equal_range(const key_type& key)             { return rep.equal_range(key); }

    std::pair<const_iterator, const_iterator>
    equal_range(const key_type& key) const       { return rep.equal_range(key); }

    mapped_type& at(const key_type& key)
    {
        iterator it = rep.find(key);
        if (it == rep.end())
            throw_exception(std::out_of_range("at: key not present"));
        return it->second;
    }

    const mapped_type& at(const key_type& key) const
    {
        const_iterator it = rep.find(key);
        if (it == rep.cend())
            throw_exception(std::out_of_range("at: key not present"));
        return it->second;
    }

#if !defined(SPP_NO_CXX11_VARIADIC_TEMPLATES)
    template <class... Args>
    std::pair<iterator, bool> emplace(Args&&... args)
    {
        return rep.emplace(std::forward<Args>(args)...);
    }

    template <class... Args>
    iterator emplace_hint(const_iterator , Args&&... args)
    {
        return rep.emplace(std::forward<Args>(args)...).first;
    }
#endif

    // Insert
    // ------
    std::pair<iterator, bool>
    insert(const value_type& obj)                    { return rep.insert(obj); }

#if !defined(SPP_NO_CXX11_RVALUE_REFERENCES)
    template< class P >
    std::pair<iterator, bool> insert(P&& obj)        { return rep.insert(std::forward<P>(obj)); }
#endif

    template <class InputIterator>
    void insert(InputIterator f, InputIterator l)    { rep.insert(f, l); }

    void insert(const_iterator f, const_iterator l)  { rep.insert(f, l); }

    iterator insert(iterator /*unused*/, const value_type& obj) { return insert(obj).first; }
    iterator insert(const_iterator /*unused*/, const value_type& obj) { return insert(obj).first; }

    // Deleted key routines - just to keep google test framework happy
    // we don't actually use the deleted key
    // ---------------------------------------------------------------
    void set_deleted_key(const key_type& key)   { rep.set_deleted_key(key); }
    void clear_deleted_key()                    { rep.clear_deleted_key();  }
    key_type deleted_key() const                { return rep.deleted_key(); }

    // Erase
    // -----
    size_type erase(const key_type& key)               { return rep.erase(key); }
    iterator  erase(iterator it)                       { return rep.erase(it); }
    iterator  erase(iterator f, iterator l)            { return rep.erase(f, l); }
    iterator  erase(const_iterator it)                 { return rep.erase(it); }
    iterator  erase(const_iterator f, const_iterator l){ return rep.erase(f, l); }

    // Comparison
    // ----------
    bool operator==(const sparse_hash_map& hs) const   { return rep == hs.rep; }
    bool operator!=(const sparse_hash_map& hs) const   { return rep != hs.rep; }


    // I/O -- this is an add-on for writing metainformation to disk
    //
    // For maximum flexibility, this does not assume a particular
    // file type (though it will probably be a FILE *).  We just pass
    // the fp through to rep.

    // If your keys and values are simple enough, you can pass this
    // serializer to serialize()/unserialize().  "Simple enough" means
    // value_type is a POD type that contains no pointers.  Note,
    // however, we don't try to normalize endianness.
    // ---------------------------------------------------------------
    typedef typename ht::NopointerSerializer NopointerSerializer;

    // serializer: a class providing operator()(OUTPUT*, const value_type&)
    //    (writing value_type to OUTPUT).  You can specify a
    //    NopointerSerializer object if appropriate (see above).
    // fp: either a FILE*, OR an ostream*/subclass_of_ostream*, OR a
    //    pointer to a class providing size_t Write(const void*, size_t),
    //    which writes a buffer into a stream (which fp presumably
    //    owns) and returns the number of bytes successfully written.
    //    Note basic_ostream<not_char> is not currently supported.
    // ---------------------------------------------------------------
    template <typename ValueSerializer, typename OUTPUT>
    bool serialize(ValueSerializer serializer, OUTPUT* fp)
    {
        return rep.serialize(serializer, fp);
    }

    // serializer: a functor providing operator()(INPUT*, value_type*)
    //    (reading from INPUT and into value_type).  You can specify a
    //    NopointerSerializer object if appropriate (see above).
    // fp: either a FILE*, OR an istream*/subclass_of_istream*, OR a
    //    pointer to a class providing size_t Read(void*, size_t),
    //    which reads into a buffer from a stream (which fp presumably
    //    owns) and returns the number of bytes successfully read.
    //    Note basic_istream<not_char> is not currently supported.
    // NOTE: Since value_type is std::pair<const Key, T>, ValueSerializer
    // may need to do a const cast in order to fill in the key.
    // NOTE: if Key or T are not POD types, the serializer MUST use
    // placement-new to initialize their values, rather than a normal
    // equals-assignment or similar.  (The value_type* passed into the
    // serializer points to garbage memory.)
    // ---------------------------------------------------------------
    template <typename ValueSerializer, typename INPUT>
    bool unserialize(ValueSerializer serializer, INPUT* fp)
    {
        return rep.unserialize(serializer, fp);
    }

    // The four methods below are DEPRECATED.
    // Use serialize() and unserialize() for new code.
    // -----------------------------------------------
    template <typename OUTPUT>
    bool write_metadata(OUTPUT *fp)       { return rep.write_metadata(fp); }

    template <typename INPUT>
    bool read_metadata(INPUT *fp)         { return rep.read_metadata(fp); }

    template <typename OUTPUT>
    bool write_nopointer_data(OUTPUT *fp) { return rep.write_nopointer_data(fp); }

    template <typename INPUT>
    bool read_nopointer_data(INPUT *fp)   { return rep.read_nopointer_data(fp); }


private:
    // The actual data
    // ---------------
    ht rep;
};

//  ----------------------------------------------------------------------
//                   S P A R S E _ H A S H _ S E T
//  ----------------------------------------------------------------------

template <class Value,
          class HashFcn  = spp_hash<Value>,
          class EqualKey = std::equal_to<Value>,
          class Alloc    = SPP_DEFAULT_ALLOCATOR<Value> >
class sparse_hash_set
{
private:
    // Apparently identity is not stl-standard, so we define our own
    struct Identity
    {
        typedef const Value& result_type;
        inline const Value& operator()(const Value& v) const { return v; }
    };

    struct SetKey
    {
        inline void operator()(Value* value, const Value& new_key) const
        {
            *value = new_key;
        }
    };

    typedef sparse_hashtable<Value, Value, HashFcn, Identity, SetKey,
                             EqualKey, Alloc> ht;

public:
    typedef typename ht::key_type              key_type;
    typedef typename ht::value_type            value_type;
    typedef typename ht::hasher                hasher;
    typedef typename ht::key_equal             key_equal;
    typedef Alloc                              allocator_type;

    typedef typename ht::size_type             size_type;
    typedef typename ht::difference_type       difference_type;
    typedef typename ht::const_pointer         pointer;
    typedef typename ht::const_pointer         const_pointer;
    typedef typename ht::const_reference       reference;
    typedef typename ht::const_reference       const_reference;

    typedef typename ht::const_iterator        iterator;
    typedef typename ht::const_iterator        const_iterator;
    typedef typename ht::const_local_iterator  local_iterator;
    typedef typename ht::const_local_iterator  const_local_iterator;


    // Iterator functions -- recall all iterators are const
    iterator       begin() const             { return rep.begin(); }
    iterator       end() const               { return rep.end(); }
    const_iterator cbegin() const            { return rep.cbegin(); }
    const_iterator cend() const              { return rep.cend(); }

    // These come from tr1's unordered_set. For us, a bucket has 0 or 1 elements.
    local_iterator begin(size_type i) const  { return rep.begin(i); }
    local_iterator end(size_type i) const    { return rep.end(i); }
    local_iterator cbegin(size_type i) const { return rep.cbegin(i); }
    local_iterator cend(size_type i) const   { return rep.cend(i); }


    // Accessor functions
    // ------------------
    allocator_type get_allocator() const     { return rep.get_allocator(); }
    hasher         hash_funct() const        { return rep.hash_funct(); }
    hasher         hash_function() const     { return hash_funct(); }  // tr1 name
    key_equal      key_eq() const            { return rep.key_eq(); }


    // Constructors
    // ------------
    explicit sparse_hash_set(size_type n = 0,
                             const hasher& hf = hasher(),
                             const key_equal& eql = key_equal(),
                             const allocator_type& alloc = allocator_type()) :
        rep(n, hf, eql, Identity(), SetKey(), alloc)
    {
    }

    explicit sparse_hash_set(const allocator_type& alloc) :
        rep(0, hasher(), key_equal(), Identity(), SetKey(), alloc)
    {
    }

    sparse_hash_set(size_type n, const allocator_type& alloc) :
        rep(n, hasher(), key_equal(), Identity(), SetKey(), alloc)
    {
    }

    sparse_hash_set(size_type n, const hasher& hf,
                    const allocator_type& alloc) :
        rep(n, hf, key_equal(), Identity(), SetKey(), alloc)
    {
    }

    template <class InputIterator>
    sparse_hash_set(InputIterator f, InputIterator l,
                    size_type n = 0,
                    const hasher& hf = hasher(),
                    const key_equal& eql = key_equal(),
                    const allocator_type& alloc = allocator_type())
        : rep(n, hf, eql, Identity(), SetKey(), alloc)
    {
        rep.insert(f, l);
    }

    template <class InputIterator>
    sparse_hash_set(InputIterator f, InputIterator l,
                    size_type n, const allocator_type& alloc)
        : rep(n, hasher(), key_equal(), Identity(), SetKey(), alloc)
    {
        rep.insert(f, l);
    }

    template <class InputIterator>
    sparse_hash_set(InputIterator f, InputIterator l,
                    size_type n, const hasher& hf, const allocator_type& alloc)
        : rep(n, hf, key_equal(), Identity(), SetKey(), alloc)
    {
        rep.insert(f, l);
    }

    sparse_hash_set(const sparse_hash_set &o) :
        rep(o.rep)
    {}

    sparse_hash_set(const sparse_hash_set &o,
                    const allocator_type& alloc) :
        rep(o.rep, alloc)
    {}

#if !defined(SPP_NO_CXX11_RVALUE_REFERENCES)
    sparse_hash_set(sparse_hash_set &&o) :
        rep(std::move(o.rep))
    {}

    sparse_hash_set(sparse_hash_set &&o,
                    const allocator_type& alloc) :
        rep(std::move(o.rep), alloc)
    {}
#endif

#if !defined(SPP_NO_CXX11_HDR_INITIALIZER_LIST)
    sparse_hash_set(std::initializer_list<value_type> init,
                    size_type n = 0,
                    const hasher& hf = hasher(),
                    const key_equal& eql = key_equal(),
                    const allocator_type& alloc = allocator_type()) :
        rep(n, hf, eql, Identity(), SetKey(), alloc)
    {
        rep.insert(init.begin(), init.end());
    }

    sparse_hash_set(std::initializer_list<value_type> init,
                    size_type n, const allocator_type& alloc) :
        rep(n, hasher(), key_equal(), Identity(), SetKey(), alloc)
    {
        rep.insert(init.begin(), init.end());
    }

    sparse_hash_set(std::initializer_list<value_type> init,
                    size_type n, const hasher& hf,
                    const allocator_type& alloc) :
        rep(n, hf, key_equal(), Identity(), SetKey(), alloc)
    {
        rep.insert(init.begin(), init.end());
    }

    sparse_hash_set& operator=(std::initializer_list<value_type> init)
    {
        rep.clear();
        rep.insert(init.begin(), init.end());
        return *this;
    }

    void insert(std::initializer_list<value_type> init)
    {
        rep.insert(init.begin(), init.end());
    }

#endif

    sparse_hash_set& operator=(const sparse_hash_set &o)
    {
        rep = o.rep;
        return *this;
    }

    void clear()                        { rep.clear(); }
    void swap(sparse_hash_set& hs)      { rep.swap(hs.rep); }


    // Functions concerning size
    // -------------------------
    size_type size() const              { return rep.size(); }
    size_type max_size() const          { return rep.max_size(); }
    bool empty() const                  { return rep.empty(); }
    size_type bucket_count() const      { return rep.bucket_count(); }
    size_type max_bucket_count() const  { return rep.max_bucket_count(); }

    size_type bucket_size(size_type i) const    { return rep.bucket_size(i); }
    size_type bucket(const key_type& key) const { return rep.bucket(key); }

    float     load_factor() const       { return size() * 1.0f / bucket_count(); }

    float max_load_factor() const      { return rep.get_enlarge_factor(); }
    void  max_load_factor(float grow)  { rep.set_enlarge_factor(grow); }

    float min_load_factor() const      { return rep.get_shrink_factor(); }
    void  min_load_factor(float shrink){ rep.set_shrink_factor(shrink); }

    void set_resizing_parameters(float shrink, float grow)
    {
        rep.set_resizing_parameters(shrink, grow);
    }

    void resize(size_type cnt)        { rep.resize(cnt); }
    void rehash(size_type cnt)        { resize(cnt); } // c++11 name
    void reserve(size_type cnt)       { resize(cnt); } // c++11

    // Lookup
    // ------
    iterator find(const key_type& key) const     { return rep.find(key); }
    bool contains(const key_type& key) const     { return rep.find(key) != rep.end(); }

    size_type count(const key_type& key) const   { return rep.count(key); }

    std::pair<iterator, iterator>
    equal_range(const key_type& key) const       { return rep.equal_range(key); }

#if !defined(SPP_NO_CXX11_VARIADIC_TEMPLATES)
    template <class... Args>
    std::pair<iterator, bool> emplace(Args&&... args)
    {
        return rep.emplace(std::forward<Args>(args)...);
    }

    template <class... Args>
    iterator emplace_hint(const_iterator , Args&&... args)
    {
        return rep.emplace(std::forward<Args>(args)...).first;
    }
#endif

    // Insert
    // ------
    std::pair<iterator, bool> insert(const value_type& obj)
    {
        std::pair<typename ht::iterator, bool> p = rep.insert(obj);
        return std::pair<iterator, bool>(p.first, p.second);   // const to non-const
    }

#if !defined(SPP_NO_CXX11_RVALUE_REFERENCES)
    template<class P>
    std::pair<iterator, bool> insert(P&& obj)        { return rep.insert(std::forward<P>(obj)); }
#endif

    template <class InputIterator>
    void insert(InputIterator f, InputIterator l)    { rep.insert(f, l); }

    void insert(const_iterator f, const_iterator l)  { rep.insert(f, l); }

    iterator insert(iterator /*unused*/, const value_type& obj) { return insert(obj).first; }

    // Deleted key - do nothing - just to keep google test framework happy
    // -------------------------------------------------------------------
    void set_deleted_key(const key_type& key) { rep.set_deleted_key(key); }
    void clear_deleted_key()                  { rep.clear_deleted_key();  }
    key_type deleted_key() const              { return rep.deleted_key(); }

    // Erase
    // -----
    size_type erase(const key_type& key)      { return rep.erase(key); }
    iterator  erase(iterator it)              { return rep.erase(it); }
    iterator  erase(iterator f, iterator l)   { return rep.erase(f, l); }

    // Comparison
    // ----------
    bool operator==(const sparse_hash_set& hs) const { return rep == hs.rep; }
    bool operator!=(const sparse_hash_set& hs) const { return rep != hs.rep; }


    // I/O -- this is an add-on for writing metainformation to disk
    //
    // For maximum flexibility, this does not assume a particular
    // file type (though it will probably be a FILE *).  We just pass
    // the fp through to rep.

    // If your keys and values are simple enough, you can pass this
    // serializer to serialize()/unserialize().  "Simple enough" means
    // value_type is a POD type that contains no pointers.  Note,
    // however, we don't try to normalize endianness.
    // ---------------------------------------------------------------
    typedef typename ht::NopointerSerializer NopointerSerializer;

    // serializer: a class providing operator()(OUTPUT*, const value_type&)
    //    (writing value_type to OUTPUT).  You can specify a
    //    NopointerSerializer object if appropriate (see above).
    // fp: either a FILE*, OR an ostream*/subclass_of_ostream*, OR a
    //    pointer to a class providing size_t Write(const void*, size_t),
    //    which writes a buffer into a stream (which fp presumably
    //    owns) and returns the number of bytes successfully written.
    //    Note basic_ostream<not_char> is not currently supported.
    // ---------------------------------------------------------------
    template <typename ValueSerializer, typename OUTPUT>
    bool serialize(ValueSerializer serializer, OUTPUT* fp)
    {
        return rep.serialize(serializer, fp);
    }

    // serializer: a functor providing operator()(INPUT*, value_type*)
    //    (reading from INPUT and into value_type).  You can specify a
    //    NopointerSerializer object if appropriate (see above).
    // fp: either a FILE*, OR an istream*/subclass_of_istream*, OR a
    //    pointer to a class providing size_t Read(void*, size_t),
    //    which reads into a buffer from a stream (which fp presumably
    //    owns) and returns the number of bytes successfully read.
    //    Note basic_istream<not_char> is not currently supported.
    // NOTE: Since value_type is const Key, ValueSerializer
    // may need to do a const cast in order to fill in the key.
    // NOTE: if Key is not a POD type, the serializer MUST use
    // placement-new to initialize its value, rather than a normal
    // equals-assignment or similar.  (The value_type* passed into
    // the serializer points to garbage memory.)
    // ---------------------------------------------------------------
    template <typename ValueSerializer, typename INPUT>
    bool unserialize(ValueSerializer serializer, INPUT* fp)
    {
        return rep.unserialize(serializer, fp);
    }

    // The four methods below are DEPRECATED.
    // Use serialize() and unserialize() for new code.
    // -----------------------------------------------
    template <typename OUTPUT>
    bool write_metadata(OUTPUT *fp)       { return rep.write_metadata(fp); }

    template <typename INPUT>
    bool read_metadata(INPUT *fp)         { return rep.read_metadata(fp); }

    template <typename OUTPUT>
    bool write_nopointer_data(OUTPUT *fp) { return rep.write_nopointer_data(fp); }

    template <typename INPUT>
    bool read_nopointer_data(INPUT *fp)   { return rep.read_nopointer_data(fp); }

private:
    // The actual data
    // ---------------
    ht rep;
};

} // spp_ namespace


// We need a global swap for all our classes as well
// -------------------------------------------------

template <class T, class Alloc>
inline void swap(spp_::sparsegroup<T,Alloc> &x, spp_::sparsegroup<T,Alloc> &y)
{
    x.swap(y);
}

template <class T, class Alloc>
inline void swap(spp_::sparsetable<T,Alloc> &x, spp_::sparsetable<T,Alloc> &y)
{
    x.swap(y);
}

template <class V, class K, class HF, class ExK, class SetK, class EqK, class A>
inline void swap(spp_::sparse_hashtable<V,K,HF,ExK,SetK,EqK,A> &x,
                 spp_::sparse_hashtable<V,K,HF,ExK,SetK,EqK,A> &y)
{
    x.swap(y);
}

template <class Key, class T, class HashFcn, class EqualKey, class Alloc>
inline void swap(spp_::sparse_hash_map<Key, T, HashFcn, EqualKey, Alloc>& hm1,
                 spp_::sparse_hash_map<Key, T, HashFcn, EqualKey, Alloc>& hm2)
{
    hm1.swap(hm2);
}

template <class Val, class HashFcn, class EqualKey, class Alloc>
inline void swap(spp_::sparse_hash_set<Val, HashFcn, EqualKey, Alloc>& hs1,
                 spp_::sparse_hash_set<Val, HashFcn, EqualKey, Alloc>& hs2)
{
    hs1.swap(hs2);
}

#endif // sparsepp_h_guard_
