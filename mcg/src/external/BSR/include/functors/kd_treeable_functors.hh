/*
 * kd-treeable functors.
 */
#ifndef FUNCTORS__KD_TREEABLE_FUNCTORS_HH
#define FUNCTORS__KD_TREEABLE_FUNCTORS_HH

#include "interfaces/kd_treeable.hh"

namespace functors {
/*
 * Imports.
 */
using interfaces::kd_tree_key;
using interfaces::kd_treeable;

/*
 * Abstract functor for kd-tree key distance operations.
 * Note that the functor is the same regardless of the type const qualifier.
 */
template <typename T, typename V = double>
class kd_treeable_key_distance_functor {
public:
   virtual ~kd_treeable_key_distance_functor() { }
   virtual V operator()(const T&, const kd_tree_key<V>&) const = 0;
};

template <typename T, typename V>
class kd_treeable_key_distance_functor<const T,V>
 : public kd_treeable_key_distance_functor<T,V> { };

/*
 * Functor for computing distance from a kd-treeable object to a kd-tree key.
 * Note that the functor is the same regardless of the type const qualifier.
 */
template <typename T, typename V = double>
class kd_tree_key_distance_functor
 : public kd_treeable_key_distance_functor<const T,V> {
public:
   V operator()(const T& t, const kd_tree_key<V>& k) const {
      return kd_treeable<T,V>::distance(t, k);
   }
};

template <typename T, typename V>
class kd_tree_key_distance_functor<const T,V>
 : public kd_tree_key_distance_functor<T,V> { };

/*
 * Globally accessible set of default kd-tree distance functors.
 */
template <typename T, typename V = double>
class kd_tree_key_distance_functors {
public:
   static const kd_tree_key_distance_functor<T,V>& f_key_distance();
};

template <typename T, typename V>
const kd_tree_key_distance_functor<T,V>& 
kd_tree_key_distance_functors<T,V>::f_key_distance() {
   static const kd_tree_key_distance_functor<T,V>* f
      = new kd_tree_key_distance_functor<T,V>();
   return *f;
}

} /* namespace functors */

#endif
