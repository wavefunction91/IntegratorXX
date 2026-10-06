// A second translation unit that includes the generator headers, so that the
// property this executable exists to check -- that a header-only build links
// from any number of translation units -- is actually exercised. In header-only
// mode these headers carry their own definitions, and any definition that loses
// its inline linkage collides with spherical_generator.cxx at link time.
//
// This has to sit in the executable rather than in integratorxx_common_ut: an
// object in a static library is pulled in only when something references it, so
// a duplicate definition there would go unnoticed.
//
// Including the headers is the whole test. A non-inline function with external
// linkage is emitted whether or not it is called, so no reference is needed.

#include <integratorxx/generators/spherical_factory.hpp>
#include <integratorxx/generators/radial_factory.hpp>
#include <integratorxx/generators/s2_factory.hpp>
