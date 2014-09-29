/*
 * Mappable interface.
 */
#ifndef INTERFACES__MAPPABLE_HH
#define INTERFACES__MAPPABLE_HH

namespace interfaces {

/*
 * Template class for mappables.
 * If a class extends mappable<T,U> then it must define a method (map_apply)
 * that can be used to map objects of type T to objects of type U.
 */
template <typename T, typename U>
class mappable {
public:
   /*
    * Destructor.
    */
   virtual ~mappable();
   
   /*
    * Map method.
    */
   virtual U& map_apply(T&) const = 0;
};

template <typename T, typename U>
mappable<T,U>::~mappable() { }

/*
 * A mappable over const T also serves as a mappable over T.
 */
template <typename T, typename U>
class mappable<const T, U> : public mappable<T,U> {
public:
   /*
    * Map method.
    * Call the map_apply version for const items.
    */
   U& map_apply(T&) const;

   /*
    * Map method.
    */
   virtual U& map_apply(const T&) const = 0;
};

template <typename T, typename U>
U& mappable<const T, U>::map_apply(T& t) const {
   const T& t_const = t;
   return (this->map_apply(t_const));
}

} /* namespace interfaces */

#endif
