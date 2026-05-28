import subprocess as sp
import os
from pathlib import Path
from . import pbs_tools as pbs

import numpy as np
from scipy.sparse import coo_matrix

import matplotlib.pyplot as plt


def run_openmx( prefix, workdir, queue='cray',
                nodes=1, ppn=1, verbose=False,
                testing=False ):

    status_file = 'scf_status.txt'

    homedir = os.getcwd()

    os.chdir( workdir )

    has_status = os.path.exists( status_file )

    node_str = f'nodes={nodes}:ppn={ppn}' 

    bash_script = """#!/bin/bash
#PBS -q batch
#PBS -N omx-{prefix}
#PBS -M enrique_blair@baylor.edu
#PBS -m abe
#PBS -l nodes=1:ppn=1

module load openmx/3.9   # enables access to openmx

# The following sets the job submission directory
#   as the working directory during job execution
cd $PBS_O_WORKDIR

# Tell mpirun the total number of MPI processes to create.                
NUMMPIPROCS=1   # Number of nodes * number of MPI processes per node (ppn$

# Specify threads per processor                                           
NUMTHREADSPERMPIPROC=1

# The following invokes OpenMX
echo Commenced job. $(date) >> {status_file}
mpirun -np $NUMMPIPROCS openmx {prefix}.dat -nt $NUMTHREADSPERMPIPROC > outputFile.std
echo Job ended. $(date) >> {status_file}""".format(**locals())

    cmdlist = ['qsub',
               '-q', queue,
               '-N', prefix, # set the job name
               '-l', node_str]

    if not has_status:
        if verbose:
            print('No status file detected. Launching calculation.')

        if not testing:

            # Form the subprocess.Popen object as a pipe
            PBS_JOB = sp.Popen(cmdlist,
                               stdin=sp.PIPE,
                               stdout=sp.PIPE,
                               stderr=sp.PIPE,
                               universal_newlines=True)

            out, err = PBS_JOB.communicate(bash_script)

            job_id = out[:7] # remove unwanted content
            pbs.start_status_file( status_file, job_id )

            print(f'\nOutput:\n' + '='*7 + f'\n{out}')

            if len(err) > 0:
                print(f'\nErrors:\n' + '='*7 + f'\n{err}')

        else:
            print('run_openmx: simulated run (testing mode).')


    else:
        if verbose:
            job_status = pbs.PBS_job_status( status_file )
            print('Status file detected. Calculation skipped.')

    os.chdir( homedir )


class OpenMXScfout:

    def __init__(self, filename, verbose=False):
        self.filename = filename
        self.read( verbose=verbose )

    def read(self, verbose=False):

        with open(self.filename, "rb") as f:

            # header
            header = np.fromfile(f, dtype=np.int32, count=6)

            if verbose:
                print(header)

            self.atomnum = header[0]
            spin_version = header[1]

            self.spinP_switch = spin_version % 4
            self.version = spin_version // 4

            self.Catomnum = header[2]
            self.Latomnum = header[3]
            self.Ratomnum = header[4]
            self.TCpyCell = header[5]

            # order_max
            self.order_max = np.fromfile(f, dtype=np.int32, count=1)[0]

            # translation vectors
            raw_atv = np.fromfile(
                f,
                dtype=np.float64,
                count=(self.TCpyCell + 1) * 4
            ).reshape(self.TCpyCell + 1, 4)

            #
            # Drop unused first column
            #
            self.atv = raw_atv[:, 1:]

            raw_atv_ijk = np.fromfile(
                f,
                dtype=np.int32,
                count=(self.TCpyCell + 1) * 4
            ).reshape(self.TCpyCell + 1, 4)

            #
            # Drop unused first column
            #
            self.atv_ijk = raw_atv_ijk[:, 1:]

            # orbitals per atom
            p_vec = np.fromfile(f, dtype=np.int32, count=self.atomnum)
            self.Total_NumOrbs = np.zeros(self.atomnum + 1, dtype=int)
            self.Total_NumOrbs[1:] = p_vec

            # neighbor counts
            p_vec = np.fromfile(f, dtype=np.int32, count=self.atomnum)
            self.FNAN = np.zeros(self.atomnum + 1, dtype=int)
            self.FNAN[1:] = p_vec

            # neighbor atom indices
            self.natn = []
            for ct_AN in range(self.atomnum + 1):
                if ct_AN == 0:
                    self.natn.append(None)
                else:
                    self.natn.append(
                        np.fromfile(
                            f, dtype=np.int32, count=self.FNAN[ct_AN] + 1
                        )
                    )

            # cell indices
            self.ncn = []
            for ct_AN in range(self.atomnum + 1):
                if ct_AN == 0:
                    self.ncn.append(None)
                else:
                    self.ncn.append(
                        np.fromfile(
                            f, dtype=np.int32, count=self.FNAN[ct_AN] + 1
                        )
                    )

            # lattice vectors
            raw_tv = np.fromfile(
                f,
                dtype=np.float64,
                count=12
            ).reshape(3, 4)

            raw_rtv = np.fromfile(
                f,
                dtype=np.float64,
                count=12
            ).reshape(3, 4)

            #
            # Drop unused first column
            #
            self.tv = raw_tv[:, 1:]
            self.rtv = raw_rtv[:, 1:]

            # atomic coordinates
            self.Gxyz = np.zeros((self.atomnum + 1, 3))

            for ct_AN in range(1, self.atomnum + 1):

                raw_xyz = np.fromfile(
                    f,
                    dtype=np.float64,
                    count=4
                )

                #
                # Drop unused first element
                #
                self.Gxyz[ct_AN] = raw_xyz[1:]

            # build orbital offsets
            self.orbital_offset = np.zeros(self.atomnum + 2, dtype=int)

            for i in range(1, self.atomnum + 1):
                self.orbital_offset[i] = (
                    self.orbital_offset[i - 1] + self.Total_NumOrbs[i - 1]
                )

            self.nbasis = self.orbital_offset[self.atomnum] + self.Total_NumOrbs[self.atomnum]

            # read the real-space Hamiltonian
            self.H_R = self.read_hamiltonian(f)

            # read the real-space overlap matrix, S
            self.S_R = self.read_overlap(f)

    def read_hamiltonian(self, f):

        #
        # Store only spin-up channel for now
        #
        rows = []
        cols = []
        data = []
        cell_indices = []

        nspin = self.spinP_switch + 1

        for spin in range(nspin):

            keep_this_spin = (spin == 0)

            for ct_AN in range(1, self.atomnum + 1):

                TNO1 = self.Total_NumOrbs[ct_AN]
                i0 = self.orbital_offset[ct_AN]

                for h_AN in range(self.FNAN[ct_AN] + 1):

                    Gh_AN = self.natn[ct_AN][h_AN]

                    Rn = self.ncn[ct_AN][h_AN]

                    TNO2 = self.Total_NumOrbs[Gh_AN]
                    j0 = self.orbital_offset[Gh_AN]

                    for i in range(TNO1):

                        vals = np.fromfile(
                            f,
                            dtype=np.float64,
                            count=TNO2
                        )

                        #
                        # Keep only spin-up
                        #
                        if keep_this_spin:

                            for j in range(TNO2):

                                rows.append(i0 + i)
                                cols.append(j0 + j)
                                data.append(vals[j])

                                cell_indices.append(Rn)

        self.H_R_rows = np.array(rows)
        self.H_R_cols = np.array(cols)
        self.H_R_cell_indices = np.array(cell_indices)
        self.H_R_data = np.array(data)

        H_R = coo_matrix(
            (data, (rows, cols)),
            shape=(self.nbasis, self.nbasis)
        )

        return H_R.tocsr()

    def read_overlap(self, f):

        rows = []
        cols = []
        data = []
        cell_indices = []

        for ct_AN in range(1, self.atomnum + 1):

            TNO1 = self.Total_NumOrbs[ct_AN]
            i0 = self.orbital_offset[ct_AN]

            for h_AN in range(self.FNAN[ct_AN] + 1):

                Gh_AN = self.natn[ct_AN][h_AN]

                Rn = self.ncn[ct_AN][h_AN]

                TNO2 = self.Total_NumOrbs[Gh_AN]
                j0 = self.orbital_offset[Gh_AN]

                for i in range(TNO1):

                    vals = np.fromfile(
                        f,
                        dtype=np.float64,
                        count=TNO2
                    )

                    for j in range(TNO2):

                        rows.append(i0 + i)
                        cols.append(j0 + j)
                        data.append(vals[j])

                        cell_indices.append(Rn)

        self.S_R_rows = np.array(rows)
        self.S_R_cols = np.array(cols)
        self.S_R_cell_indices = np.array(cell_indices)
        self.S_R_data = np.array(data)

        S_R = coo_matrix(
            (data, (rows, cols)),
            shape=(self.nbasis, self.nbasis)
        )

        return S_R.tocsr()



    def Hk(self, kvec):

        """
        Construct Bloch Hamiltonian H(k).

        Parameters
        ----------
        kvec : array-like shape (3,)
        k-vector in Cartesian coordinates (Bohr^-1)

        Returns
        -------
        Hk : ndarray (nbasis, nbasis)
        Complex Bloch Hamiltonian
        """

        Hk = np.zeros(
            (self.nbasis, self.nbasis),
            dtype=np.complex128
        )

        for row, col, cell_index, val in zip(
                self.H_R_rows,
                self.H_R_cols,
                self.H_R_cell_indices,
                self.H_R_data
        ):

            #
            # Real-space lattice translation vector
            #
            R = self.atv[cell_index]

            #
            # Bloch phase
            #
            phase = np.exp(1j * np.dot(kvec, R))

            Hk[row, col] += val * phase

        return Hk

    def Sk(self, kvec):

        """
        Construct Bloch overlap matrix S(k).
        """

        Sk = np.zeros(
            (self.nbasis, self.nbasis),
            dtype=np.complex128
        )

        for row, col, cell_index, val in zip(
                self.S_R_rows,
                self.S_R_cols,
                self.S_R_cell_indices,
                self.S_R_data
        ):

            R = self.atv[cell_index]

            phase = np.exp(
                1j * np.dot(kvec, R)
            )

            Sk[row, col] += val * phase

        return Sk

    def reciprocal_vectors(self):
        """
        Reciprocal lattice vectors in Cartesian coordinates.
        """

        a1 = self.tv[0]
        a2 = self.tv[1]
        a3 = self.tv[2]

        volume = np.dot(a1, np.cross(a2, a3))

        b1 = 2*np.pi * np.cross(a2, a3) / volume
        b2 = 2*np.pi * np.cross(a3, a1) / volume
        b3 = 2*np.pi * np.cross(a1, a2) / volume

        return np.array([b1, b2, b3])

def read_openmx_band(filename, reciprocal_vectors):
    """
    Read OpenMX .Band file.

    Returns
    -------
    kdist : ndarray
        1D k coordinate

    bands : ndarray
        shape (nk, nbands), in eV

    fermi_energy : float
        Fermi energy in eV
    """

    with open(filename, "r") as f:

        lines = [
            line.strip()
            for line in f
            if line.strip()
        ]

    #
    # First line:
    #
    # nbands nspin Ef
    #
    tokens = lines[0].split()

    nbands = int(tokens[0])

    nspin = int(tokens[1])

    Ef = float(tokens[2]) * 27.2114

    #
    # Number of path segments
    #
    nkpath = int(lines[2])

    #
    # Skip path definitions
    #
    line_index = 3 + nkpath

    kdist = []

    bands = []

    cumulative_k = 0.0

    prev_kvec = None

    while line_index < len(lines):

        #
        # Header line:
        #
        # nbands kx ky kz
        #
        header = lines[line_index].split()

        if len(header) != 4:
            break

        nb = int(header[0])

        kx = float(header[1])
        ky = float(header[2])
        kz = float(header[3])

        kfrac = np.array([kx, ky, kz])

        kvec = (
            kfrac[0] * reciprocal_vectors[0]
            + kfrac[1] * reciprocal_vectors[1]
            + kfrac[2] * reciprocal_vectors[2]
        )

        #
        # Eigenvalue line
        #
        eigvals = np.array(
            [float(x) for x in lines[line_index + 1].split()]
        )

        #
        # Hartree -> eV
        #
        eigvals *= 27.2114

        #
        # Build cumulative k-distance
        #
        if prev_kvec is not None:

            cumulative_k += np.linalg.norm(
                kvec - prev_kvec
            )

        kdist.append(cumulative_k)

        bands.append(eigvals)

        prev_kvec = kvec

        line_index += 2

    return (
        np.array(kdist),
        np.array(bands),
        Ef
    )
