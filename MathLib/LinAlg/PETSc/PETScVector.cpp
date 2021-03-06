/**
 * \file  PETScVector.cpp
 * \brief Definition of member functions of class PETScVector,
 *        which provides an interface to PETSc vector routines.
 *
 *   Note: the return message of PETSc routines is ommited in
 *         the source code. If it is really needed, it can be activated by
 *         adding a PetscErrorCode type variable before each PETSc fucntion
 *
 * \author Wenqing Wang
 * \date Nov 2011 - Sep 2013
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#include "PETScVector.h"

#include <algorithm>
#include <cassert>

namespace MathLib
{
PETScVector::PETScVector(const PetscInt vec_size, const bool is_global_size)
{
    _v.reset(new PETSc_Vec);

    if( is_global_size ) {
        VecCreate(PETSC_COMM_WORLD, _v.get());
        VecSetSizes(*_v, PETSC_DECIDE, vec_size);
    } else {
        // Fix size partitioning
        // the size can be associated to specific memory allocation of a matrix
        VecCreateMPI(PETSC_COMM_WORLD, vec_size, PETSC_DECIDE, _v.get());
    }

    config();
}

PETScVector::PETScVector(const PetscInt vec_size,
                         const std::vector<PetscInt>& ghost_ids,
                         const bool is_global_size)
    : _size_ghosts{ghost_ids.size()}
    , _has_ghost_id{true}
{
    _v.reset(new PETSc_Vec);

    PetscInt nghosts = static_cast<PetscInt>( ghost_ids.size() );
    if ( is_global_size )
    {
        VecCreateGhost(PETSC_COMM_WORLD, PETSC_DECIDE, vec_size, nghosts,
                       ghost_ids.data(), _v.get());
    }
    else
    {
        VecCreate(PETSC_COMM_WORLD, _v.get());
        VecSetType(*_v, VECMPI);
        VecSetSizes(*_v, vec_size, PETSC_DECIDE);
        VecMPISetGhost(*_v, nghosts, ghost_ids.data());
    }

    config();
}

PETScVector::PETScVector(const PETScVector &existing_vec, const bool deep_copy)
{
    shallowCopy(existing_vec);

    // Copy values
    if(deep_copy)
    {
        VecCopy(*existing_vec._v, *_v);
    }
}

PETScVector::PETScVector(PETScVector &&other)
    : _v{std::move(other._v)}
    , _v_loc{std::move(other._v_loc)}
    , _start_rank{other._start_rank}
    , _end_rank{other._end_rank}
    , _size{other._size}
    , _size_loc{other._size_loc}
    , _size_ghosts{other._size_ghosts}
    , _has_ghost_id{other._has_ghost_id}
{}

void PETScVector::config()
{
    VecSetFromOptions(*_v);
    // VecSetUp(*_v); // for petsc ver.>3.3
    VecGetOwnershipRange(*_v, &_start_rank, &_end_rank);

    VecGetLocalSize(*_v, &_size_loc);
    VecGetSize(*_v, &_size);

    VecSetOption(*_v, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
}

void PETScVector::finalizeAssembly()
{
    VecAssemblyBegin(*_v);
    VecAssemblyEnd(*_v);
}

void PETScVector::gatherLocalVectors( PetscScalar local_array[],
                                      PetscScalar global_array[])
{
    // Collect vectors from processors.
    int size_rank;
    MPI_Comm_size(PETSC_COMM_WORLD, &size_rank);

    // number of elements to be sent for each rank
    std::vector<PetscInt>  i_cnt(size_rank);

    MPI_Allgather(&_size_loc, 1, MPI_INT, &i_cnt[0], 1, MPI_INT, PETSC_COMM_WORLD);

    // collect local array
    PetscInt offset = 0;
    // offset in the receive vector of the data from each rank
    std::vector<PetscInt>  i_disp(size_rank);
    for(PetscInt i=0; i<size_rank; i++)
    {
        i_disp[i] = offset;
        offset += i_cnt[i];
    }

    MPI_Allgatherv(local_array, _size_loc, MPI_DOUBLE,
                   global_array, &i_cnt[0], &i_disp[0], MPI_DOUBLE, PETSC_COMM_WORLD);

}

void PETScVector::getGlobalVector(PetscScalar u[])
{

#ifdef TEST_MEM_PETSC
    PetscLogDouble mem1, mem2;
    PetscMemoryGetCurrentUsage(&mem1);
#endif

    PetscScalar *xp = nullptr;
    VecGetArray(*_v, &xp);

    gatherLocalVectors(xp, u);

    //This following line may be needed late on
    //  for a communication load balance:
    //MPI_Barrier(PETSC_COMM_WORLD);

    VecRestoreArray(*_v, &xp);

    //TEST
#ifdef TEST_MEM_PETSC
    PetscMemoryGetCurrentUsage(&mem2);
    PetscPrintf(PETSC_COMM_WORLD, "### Memory usage by Updating. Before :%f After:%f Increase:%d\n", mem1, mem2, (int)(mem2 - mem1));
#endif
}

void PETScVector::copyValues(std::vector<double>& u) const
{
    assert(u.size() == (std::size_t) (getLocalSize() + getGhostSize()));

    double* loc_x = getLocalVector();
    std::copy_n(loc_x, getLocalSize() + getGhostSize(), u.begin());
    restoreArray(loc_x);
}

PetscScalar* PETScVector::getLocalVector() const
{
    PetscScalar *loc_array;
    if (_has_ghost_id)
    {
        VecGhostUpdateBegin(*_v, INSERT_VALUES, SCATTER_FORWARD);
        VecGhostUpdateEnd(*_v, INSERT_VALUES, SCATTER_FORWARD);
        VecGhostGetLocalForm(*_v, &_v_loc);
        VecGetArray(_v_loc, &loc_array);
    }
    else
       VecGetArray(*_v, &loc_array);
    return loc_array;
}

void PETScVector::restoreArray(PetscScalar* array) const
{
    if (_has_ghost_id)
    {
        VecRestoreArray(_v_loc, &array);
        VecGhostRestoreLocalForm(*_v, &_v_loc);
    }
    else
        VecRestoreArray(*_v, &array);
}

PetscScalar PETScVector::getNorm(MathLib::VecNormType nmtype) const
{
    NormType petsc_norm = NORM_1;
    switch(nmtype)
    {
        case MathLib::VecNormType::NORM1:
            petsc_norm = NORM_1;
            break;
        case MathLib::VecNormType::NORM2:
            petsc_norm = NORM_2;
            break;
        case MathLib::VecNormType::INFINITY_N:
            petsc_norm = NORM_INFINITY;
            break;
        default:
            break;
    }

    PetscScalar norm = 0.;
    VecNorm(*_v, petsc_norm, &norm);
    return norm;
}

void PETScVector::viewer(const std::string &file_name, const PetscViewerFormat vw_format) const
{
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
    PetscViewerPushFormat(viewer, vw_format);

    PetscObjectSetName((PetscObject)*_v, file_name.c_str());
    VecView(*_v, viewer);

#define  nEXIT_TEST
#ifdef EXIT_TEST
    VecDestroy(_v);
    PetscFinalize();
    exit(0);
#endif

}

void PETScVector::shallowCopy(const PETScVector &v)
{
    destroy();

    _v.reset(new PETSc_Vec);

    VecDuplicate(*v._v, _v.get());

    _start_rank   = v._start_rank;
    _end_rank     = v._end_rank;
    _size         = v._size;
    _size_loc     = v._size_loc;
    _size_ghosts  = v._size_ghosts;
    _has_ghost_id = v._has_ghost_id;

    VecSetOption(*_v, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
}

void finalizeVectorAssembly(PETScVector &vec)
{
    vec.finalizeAssembly();
}

} //end of namespace
