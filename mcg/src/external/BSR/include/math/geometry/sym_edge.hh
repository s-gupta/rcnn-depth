/*
 * Sym-edge.
 *
 * Sym-edges are quad-edges restricted to the primary oriented graph.  Their
 * implementation does not maintain the data structures for accessing the
 * rotated, flipped, or dual graphs.  Only the symmetric counterpart of each
 * edge is stored.
 *
 * Each vertex has an associated data item and each directed edge may
 * optionally have an associated data item.
 */
#ifndef MATH__GEOMETRY__SYM_EDGE_HH
#define MATH__GEOMETRY__SYM_EDGE_HH

#include "lang/exceptions/ex_not_found.hh"
#include "lang/pointers/auto_ptr.hh"
#include "math/geometry/point_2D.hh"

namespace math {
namespace geometry {
/*
 * Imports.
 */
using lang::exceptions::ex_not_found;
using lang::pointers::auto_ptr;

/*
 * Sym-edge.
 * V is the type of data associated with each vertex.
 * E is the (optional) type of data associated with each directed edge.
 */
template <typename V = point_2D, typename E = void*>
class sym_edge {
public:
   /*
    * Constructors.
    * Optionally specify (undirected or directed) edge data.
    */
   explicit sym_edge(V& /* origin */, V& /* destination */);
   explicit sym_edge(V& /* origin */, V& /* destination */, E&     /* data */);
   explicit sym_edge(V& /* origin */, V& /* destination */, E&, E& /* data */);
   
   /*
    * Destructor.
    * Disconnect the edge from the edge structure and delete it.
    */
   virtual ~sym_edge();

   /*
    * Get/set origin.
    */
   V& origin() const;
   void origin(V&);
   
   /*
    * Get/set destination.
    */
   V& destination() const;
   void destination(V&);

   /*
    * Check whether the edge has associated data.
    */
   bool has_data() const;
   
   /*
    * Get edge data.
    * Return the data associated with the edge.
    * Throw an ex_not_found exception if there is no associated data.
    */
   E& data() const;

   /*
    * Set edge data.
    * Associate the given data item with the edge.
    */
   void data(E&);

   /*
    * Return symmetric edge (with same orientation and opposite direction).
    */
   sym_edge<V,E>& sym();
   const sym_edge<V,E>& sym() const;

   /*
    * Return the previous edge with the same origin (move clockwise).
    */
   sym_edge<V,E>& origin_prev();
   const sym_edge<V,E>& origin_prev() const;
   
   /*
    * Return the next edge with the same origin (move counterclockwise).
    */
   sym_edge<V,E>& origin_next();
   const sym_edge<V,E>& origin_next() const;

   /*
    * Return the previous edge with the same destination.
    */
   sym_edge<V,E>& dest_prev();
   const sym_edge<V,E>& dest_prev() const;
   
   /*
    * Return the next edge with the same destination.
    */
   sym_edge<V,E>& dest_next();
   const sym_edge<V,E>& dest_next() const;

   /*
    * Return the previous edge on the left face.
    */
   sym_edge<V,E>& left_prev();
   const sym_edge<V,E>& left_prev() const;

   /*
    * Return the next edge on the left face.
    */
   sym_edge<V,E>& left_next();
   const sym_edge<V,E>& left_next() const;
   
   /*
    * Return the previous edge on the right face.
    */
   sym_edge<V,E>& right_prev();
   const sym_edge<V,E>& right_prev() const;
   
   /*
    * Return the next edge on the right face.
    */
   sym_edge<V,E>& right_next();
   const sym_edge<V,E>& right_next() const;

   /*
    * Swap the edge for the diagonal edge formed by sliding its endpoints 
    * along their immediately previous edges.  Return a reference to the edge.
    */
   sym_edge<V,E>& swap();

   /*
    * Disconnect the edge from the edge structure.
    * Return a reference to the edge.
    */
   sym_edge<V,E>& disconnect();

   /*
    * Reconnect an edge that has been disconnected.  Set it to connect the two
    * given edges (from the destination of the first edge to the origin of the
    * second).  Return a reference to the edge.
    */
   sym_edge<V,E>& reconnect(sym_edge<V,E>&, sym_edge<V,E>&);
    
   /*
    * Create an edge connecting the two given edges (from the destination of
    * the first edge to the origin of the second).  Return the new edge.
    */
   static auto_ptr< sym_edge<V,E> > connect(
      sym_edge<V,E>&, sym_edge<V,E>&
   );
   
   /*
    * Create an edge connecting the two given edges and set its edge data.
    * Return the new edge.
    */
   static auto_ptr< sym_edge<V,E> > connect(
      sym_edge<V,E>&, sym_edge<V,E>&, E&
   );

   /*
    * Create an edge connecting the two given edges and set its edge data.
    * Return the new edge.
    */
   static auto_ptr< sym_edge<V,E> > connect(
      sym_edge<V,E>&, sym_edge<V,E>&, E&, E&
   );
   
   /*
    * Splice two edges.
    */
   static void splice(sym_edge<V,E>&, sym_edge<V,E>&);

protected:
   /*
    * Protected constructor.
    * Create and return an empty sym-edge.
    */
   sym_edge();

   /*
    * Protected copy constructor.
    * Sym-edges should only be copied internally.
    */
   sym_edge(const sym_edge<V,E>&);

   /*
    * Sym-edge data.
    */
   V*             _v_data;    /* vertex data */
   E*             _e_data;    /* edge data (optional) */
   sym_edge<V,E>* _sym;       /* same edge with opposite direction */
   sym_edge<V,E>* _oprev;     /* previous edge with the same origin */
   sym_edge<V,E>* _onext;     /* next edge with the same origin */
};

/***************************************************************************
 * Constructors and destructor.
 ***************************************************************************/

/*
 * Protected constructor.
 * Create and return an empty sym-edge.
 */
template <typename V, typename E>
sym_edge<V,E>::sym_edge()
 : _v_data(NULL), _e_data(NULL), _sym(NULL), _oprev(this), _onext(this)
{ }

/*
 * Constructor.
 * Create a sym-edge with the given origin and destination.
 */
template <typename V, typename E>
sym_edge<V,E>::sym_edge(V& v_origin, V& v_dest)
 : _v_data(&v_origin), _e_data(NULL), _sym(NULL), _oprev(this), _onext(this)
{
   _sym = new sym_edge<V,E>();
   _sym->_v_data = &v_dest;
   _sym->_sym = this;
}

/*
 * Constructor.
 * Create a sym-edge with the given origin and destination.
 * In addition, associate the given data item with the undirected edge.
 */
template <typename V, typename E>
sym_edge<V,E>::sym_edge(V& v_origin, V& v_dest, E& e_data)
 : _v_data(&v_origin), _e_data(&e_data), _sym(NULL), _oprev(this), _onext(this)
{
   _sym = new sym_edge<V,E>();
   _sym->_v_data = &v_dest;
   _sym->_e_data = &e_data;
   _sym->_sym = this;
}

/*
 * Constructor.
 * Create a sym-edge with the given origin and destination.
 * In addition, associate the given data items with the directed edges.
 */
template <typename V, typename E>
sym_edge<V,E>::sym_edge(V& v_origin, V& v_dest, E& e_data, E& e_sym_data)
 : _v_data(&v_origin), _e_data(&e_data), _sym(NULL), _oprev(this), _onext(this)
{
   _sym = new sym_edge<V,E>();
   _sym->_v_data = &v_dest;
   _sym->_e_data = &e_sym_data;
   _sym->_sym = this;
}

/*
 * Protected copy constructor.
 * Sym-edges should only be copied internally.
 */
template <typename V, typename E>
sym_edge<V,E>::sym_edge(const sym_edge<V,E>& e)
 : _v_data(e._vdata),
   _e_data(e._e_data),
   _sym(e._sym),
   _oprev(e._oprev),
   _onext(e._onext)
{
   /* unlink original edge */
   e._sym = NULL;
   e._oprev = &e;
   e._onext = &e;
   /* link symmetric edge */
   _sym->_sym = this;
   /* check previous and next edges for self-reference */
   if (_oprev == &e) { _oprev = this; }
   if (_onext == &e) { _onext = this; }
   /* link previous and next edges */
   _oprev->_onext = this;
   _onext->_oprev = this;
}

/*
 * Destructor.
 * Disconnect the edge from the edge structure and delete it.
 */
template <typename V, typename E>
sym_edge<V,E>::~sym_edge() {
   /* disconnect edge form edge structure */
   sym_edge<V,E>::splice(*this, *_oprev);
   /* delete symmetric edge */
   if (_sym != NULL) {
      _sym->_sym = NULL;
      delete _sym;
   }
}

/***************************************************************************
 * Vertex and edge data.
 ***************************************************************************/

/*
 * Get origin.
 */
template <typename V, typename E>
V& sym_edge<V,E>::origin() const {
   return *_v_data;
}

/*
 * Set origin.
 */
template <typename V, typename E>
void sym_edge<V,E>::origin(V& v) {
   _v_data = &v;
}

/*
 * Get destination.
 */
template <typename V, typename E>
V& sym_edge<V,E>::destination() const {
   return *(_sym->_v_data);
}

/*
 * Set destination.
 */
template <typename V, typename E>
void sym_edge<V,E>::destination(V& v) {
   _sym->_v_data = &v;
}

/*
 * Check whether the edge has associated data.
 */
template <typename V, typename E>
bool sym_edge<V,E>::has_data() const {
   return (_e_data != NULL);
}

/*
 * Get edge data.
 * Return the data associated with the edge.
 * Throw an ex_not_found exception if there is no associated data.
 */
template <typename V, typename E>
E& sym_edge<V,E>::data() const {
   if (_e_data == NULL)
      throw ex_not_found("no data associated with sym_edge");
   else
      return *_e_data;
}

/*
 * Set edge data.
 * Associate the given data item with the edge.
 */
template <typename V, typename E>
void sym_edge<V,E>::data(E& e) {
   _e_data = &e;
}

/***************************************************************************
 * Edge counterparts.
 ***************************************************************************/

/*
 * Return edge with same orientation and opposite direction.
 */
template <typename V, typename E>
sym_edge<V,E>& sym_edge<V,E>::sym() {
   return *_sym;
}

template <typename V, typename E>
const sym_edge<V,E>& sym_edge<V,E>::sym() const {
   return *_sym;
}

/*
 * Return the previous edge with the same origin (move clockwise).
 */
template <typename V, typename E>
sym_edge<V,E>& sym_edge<V,E>::origin_prev() {
   return *_oprev;
}

template <typename V, typename E>
const sym_edge<V,E>& sym_edge<V,E>::origin_prev() const {
   return *_oprev;
}

/*
 * Return the next edge with the same origin (move counterclockwise).
 */
template <typename V, typename E>
sym_edge<V,E>& sym_edge<V,E>::origin_next() {
   return *_onext;
}

template <typename V, typename E>
const sym_edge<V,E>& sym_edge<V,E>::origin_next() const {
   return *_onext;
}

/*
 * Return the previous edge with the same destination.
 */
template <typename V, typename E>
sym_edge<V,E>& sym_edge<V,E>::dest_prev() {
   return *(_sym->_oprev->_sym);
}

template <typename V, typename E>
const sym_edge<V,E>& sym_edge<V,E>::dest_prev() const {
   return *(_sym->_oprev->_sym);
}

/*
 * Return the next edge with the same destination.
 */
template <typename V, typename E>
sym_edge<V,E>& sym_edge<V,E>::dest_next() {
   return *(_sym->_onext->_sym);
}  

template <typename V, typename E>
const sym_edge<V,E>& sym_edge<V,E>::dest_next() const {
   return *(_sym->_onext->_sym);
}

/*
 * Return the previous edge on the left face.
 */
template <typename V, typename E>
sym_edge<V,E>& sym_edge<V,E>::left_prev() {
   return *(_onext->_sym);
}

template <typename V, typename E>
const sym_edge<V,E>& sym_edge<V,E>::left_prev() const {
   return *(_onext->_sym);
}

/*
 * Return the next edge on the left face.
 */
template <typename V, typename E>
sym_edge<V,E>& sym_edge<V,E>::left_next() {
   return *(_sym->_oprev);
}

template <typename V, typename E>
const sym_edge<V,E>& sym_edge<V,E>::left_next() const {
   return *(_sym->_oprev);
}

/*
 * Return the previous edge on the right face.
 */
template <typename V, typename E>
sym_edge<V,E>& sym_edge<V,E>::right_prev() {
   return *(_sym->_onext);
}

template <typename V, typename E>
const sym_edge<V,E>& sym_edge<V,E>::right_prev() const {
   return *(_sym->_onext);
}

/*
 * Return the next edge on the right face.
 */
template <typename V, typename E>
sym_edge<V,E>& sym_edge<V,E>::right_next() {
   return *(_oprev->_sym);
}
   
template <typename V, typename E>
const sym_edge<V,E>& sym_edge<V,E>::right_next() const {
   return *(_oprev->_sym);
}

/***************************************************************************
 * Edge operations.
 ***************************************************************************/

/*
 * Swap the edge for the diagonal edge formed by sliding its endpoints 
 * along their immediately previous edges.  Return a reference to the edge.
 */
template <typename V, typename E>
sym_edge<V,E>& sym_edge<V,E>::swap() {
   sym_edge<V,E>& a(*(_oprev));
   sym_edge<V,E>& b(*(_sym->_oprev));
   sym_edge<V,E>::splice(*this, a);
   sym_edge<V,E>::splice(*_sym, b);
   sym_edge<V,E>::splice(*this, a.left_next());
   sym_edge<V,E>::splice(*_sym, b.left_next());
   _v_data       = &(a.destination());
   _sym->_v_data = &(b.destination());
   return *this;
}

/*
 * Disconnect the edge from the edge structure.
 * Return a reference to the edge.
 */
template <typename V, typename E>
sym_edge<V,E>& sym_edge<V,E>::disconnect() {
   sym_edge<V,E>::splice(*this, *_oprev);
   sym_edge<V,E>::splice(*_sym, *(_sym->_oprev));
   return *this;
}

/*
 * Reconnect an edge that has been disconnected.  Set it to connect the two
 * given edges (from the destination of the first edge to the origin of the
 * second).  Return a reference to the edge.
 */
template <typename V, typename E>
sym_edge<V,E>& sym_edge<V,E>::reconnect(sym_edge<V,E>& a, sym_edge<V,E>& b) {
   /* set edge origin and destination */
   this->origin(a.destination());
   this->destination(b.origin());
   /* set edge to connect a to b */
   sym_edge<V,E>::splice(*this, a.left_next());
   sym_edge<V,E>::splice(this->sym(), b);
   return *this;
}

/*
 * Create an edge connecting the two given edges (from the destination of
 * the first edge to the origin of the second).  Return the new edge.
 */
template <typename V, typename E>
auto_ptr< sym_edge<V,E> > sym_edge<V,E>::connect(
   sym_edge<V,E>& a, sym_edge<V,E>& b)
{
   auto_ptr< sym_edge<V,E> > e(new sym_edge<V,E>(a.destination(), b.origin()));
   sym_edge<V,E>::splice(*e, a.left_next());
   sym_edge<V,E>::splice(e->sym(), b);
   return e;
}

/*
 * Create an edge connecting the two given edges and set its edge data.
 * Return the new edge.
 */
template <typename V, typename E>
auto_ptr< sym_edge<V,E> > sym_edge<V,E>::connect(
   sym_edge<V,E>& a, sym_edge<V,E>& b, E& e_data)
{
   auto_ptr< sym_edge<V,E> > e = sym_edge<V,E>::connect(a, b);
   e->data(e_data);
   return e;
}

/*
 * Create an edge connecting the two given edges and set its edge data.
 * Return the new edge.
 */
template <typename V, typename E>
auto_ptr< sym_edge<V,E> > sym_edge<V,E>::connect(
   sym_edge<V,E>& a, sym_edge<V,E>& b, E& e_data, E& e_sym_data)
{
   auto_ptr< sym_edge<V,E> > e = sym_edge<V,E>::connect(a, b);
   e->data(e_data);
   e->_sym->data(e_sym_data);
   return e;
}

/*
 * Splice two edges.
 */
template <typename V, typename E>
void sym_edge<V,E>::splice(sym_edge<V,E>& a, sym_edge<V,E>& b) {
   sym_edge<V,E>* temp = a._onext;
   a._onext = b._onext;
   b._onext->_oprev = &a;
   b._onext = temp;
   temp->_oprev = &b;
}

} /* namespace geometry */
} /* namespace math */

#endif
