/*
 * Two-column lookup table, shared by the material and the source plugin.
 *
 * Both plugins used to carry one copy of this loader per data file -- ten in
 * matlib.cpp and one in current.cpp, sixty-eight lines each, differing only in
 * the file name. They also located their data with the __FILE__ trick, which
 * baked the SOURCE directory into the binary and broke the moment the sources
 * moved into ./src and the data into ./data. Both problems are fixed here:
 * one implementation, and the data directory comes from CMake as
 * BELFEM_USER_DATA_DIR, so the library does not depend on where it is launched
 * from either.
 */

#ifndef TAPE_QUENCH_USER_TABLE_HPP
#define TAPE_QUENCH_USER_TABLE_HPP

#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <belfem_user_api>

#ifndef BELFEM_USER_DATA_DIR
#error "BELFEM_USER_DATA_DIR is not set -- see the target_compile_definitions in CMakeLists.txt"
#endif

namespace usermat
{
//------------------------------------------------------------------------------

    /**
     * @brief A two-column ( x, y ) text table with linear interpolation.
     *
     * The file is read on the first call and kept for the life of the process.
     * Arguments outside the sampled range are clamped to the end values, which
     * is what the original per-file loaders did.
     *
     * Declare one as a function-local static:
     *
     *     real cu_cp( const Material* mat, real T )
     *     {
     *         static usermat::Table tTable( "cp_Cu.txt" ) ;
     *         return tTable( T ) ;
     *     }
     */
    class Table
    {
        const std::string             mFile ;
        std::vector< belfem::real >   mX ;
        std::vector< belfem::real >   mY ;
        bool                          mLoaded = false ;

//------------------------------------------------------------------------------
    public:
//------------------------------------------------------------------------------

        explicit Table( const std::string & File ) : mFile( File )
        {
        }

//------------------------------------------------------------------------------

        belfem::real
        operator()( const belfem::real x )
        {
            if ( ! mLoaded )
            {
                this->load() ;
            }

            if ( x <= mX.front() ) return mY.front() ;
            if ( x >= mX.back()  ) return mY.back() ;

            // index of the first sample strictly greater than x. The two
            // clamps above guarantee 1 <= k <= mX.size() - 1, so neither
            // k - 1 nor k can run off the ends.
            const std::size_t k = std::distance(
                    mX.begin(), std::upper_bound( mX.begin(), mX.end(), x ) ) ;

            const belfem::real t = ( x - mX[ k-1 ] ) / ( mX[ k ] - mX[ k-1 ] ) ;

            return mY[ k-1 ] + t * ( mY[ k ] - mY[ k-1 ] ) ;
        }

//------------------------------------------------------------------------------
    private:
//------------------------------------------------------------------------------

        void
        load()
        {
            const std::string tPath = std::string( BELFEM_USER_DATA_DIR ) + mFile ;

            std::ifstream tFile( tPath ) ;

            BELFEM_ERROR( tFile.is_open(),
                "Cannot open user data file %s", tPath.c_str() ) ;

            std::string  tLine ;
            belfem::real x ;
            belfem::real y ;

            while ( std::getline( tFile, tLine ) )
            {
                std::istringstream tStream( tLine ) ;

                if ( tStream >> x >> y )
                {
                    mX.push_back( x ) ;
                    mY.push_back( y ) ;
                }
            }

            tFile.close() ;

            BELFEM_ERROR( mX.size() > 1,
                "User data file %s contains fewer than two samples", tPath.c_str() ) ;

            mLoaded = true ;
        }
    };

//------------------------------------------------------------------------------
}
#endif //TAPE_QUENCH_USER_TABLE_HPP
