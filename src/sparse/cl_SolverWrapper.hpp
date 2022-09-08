//
// Created by Christian Messe on 10.07.20.
//

#ifndef BELFEM_CL_SOLVERWRAPPER_HPP
#define BELFEM_CL_SOLVERWRAPPER_HPP

#include "typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "cl_SpMatrix.hpp"
#include "en_SolverEnums.hpp"

namespace belfem
{
    namespace solver
    {
        /**
         * parent class for solver specific data
         */
        class Wrapper
        {

//------------------------------------------------------------------------------

            // rank of this proc
            const proc_t mMyRank ;

            const proc_t mCommSize ;

            // label of this solver type
            const string mLabel ;

            // flag telling if we have been initialized
            bool mIsInitialized = false ;

//------------------------------------------------------------------------------
        public:
//------------------------------------------------------------------------------

            Wrapper( const string & aLabel );

//------------------------------------------------------------------------------

            virtual ~Wrapper() = default;

//------------------------------------------------------------------------------

            /**
             * returns the name of the solver as string
             * @return
             */
            const string &
            label() const ;

//------------------------------------------------------------------------------

            virtual void
            solve( SpMatrix & aMatrix,
                   Vector <real> & aLHS,
                   Vector <real> & aRHS );

//------------------------------------------------------------------------------

            virtual void
            solve( SpMatrix & aMatrix,
                   Matrix <real> & aLHS,
                   Matrix <real> & aRHS );

//------------------------------------------------------------------------------

            /**
             * tells if the initialitzation routine has already been called
             */
            bool
            is_initialized() const ;

//------------------------------------------------------------------------------

            /**
             * returns the communication rank of this proc
             */
            proc_t
            rank() const ;

//------------------------------------------------------------------------------

            /**
             * returns the communication size of this proc
             */
            proc_t
            comm_size() const ;

//------------------------------------------------------------------------------

        virtual void
        initialize(
                SpMatrix & aMatrix,
                const SymmetryMode aSymmetryMode = SymmetryMode::GeneralSymmetric,
                const int aNumRhsColumns = 1 );

//------------------------------------------------------------------------------

        virtual void
        free();

//------------------------------------------------------------------------------

        /**
         * returns the determinant, if supported by the solver
         * and computation was requested
         */
        virtual real
        get_determinant() const ;

//------------------------------------------------------------------------------

        /**
        * returns the conditioning numbner, if supported by the solver
        * and computation was requested
        */
        virtual real
        get_cond0() const ;

//------------------------------------------------------------------------------

        /**
         * returns the conditioning numbner, if supported by the solver
         * and computation was requested
        */
       virtual real
       get_cond1() const ;

//------------------------------------------------------------------------------
        protected:
//------------------------------------------------------------------------------

            virtual void
            initialize();

//------------------------------------------------------------------------------

            void
            mat2vec( const Matrix< real > & aM,
                           Vector< real > & aV );

//------------------------------------------------------------------------------

            void
            vec2mat( const Vector< real > & aV,
                           Matrix< real > & aM );



//------------------------------------------------------------------------------
        };

//------------------------------------------------------------------------------

        inline bool
        Wrapper::is_initialized() const
        {
            return mIsInitialized ;
        }

//------------------------------------------------------------------------------

        inline proc_t
        Wrapper::rank() const
        {
            return mMyRank ;
        }

//------------------------------------------------------------------------------

        inline proc_t
        Wrapper::comm_size() const
        {
            return mCommSize ;
        }

//------------------------------------------------------------------------------
    }
}

#endif //BELFEM_CL_SOLVERWRAPPER_HPP
