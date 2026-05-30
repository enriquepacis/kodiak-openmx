# README

## Introduction

This repository is focused on running OpenMX calculations and extracting the
real-space Hamiltonian.

## Basis Functions

### Pseudo-atomic Orbitals (PAOs)

OpenMX uses localized atomic-like basis functions called pseudo-atomic orbitals
(PAOs). Each basis function is $\phi_{i\alpha}(\mathbf{r})$. It is a numerical
atomic orbital centered on atom $i$. The angular part of $\phi_{i\alpha}$ is
defined by spherical harmonics, and the radial part is the numerical solution of
the atomic problem.

The atomic problem means solving the Schrodinger equation for an isolated atom. That is, for an atom with nuclear charge $Z$, we solve the spherically symmetric electronic structure problem,
$$\hat{H}_{\mbox{atom}} \psi_{n\ell m} = E_{n\ell} \psi_{n\ell m}$$
OpenMX does not use analytic basis functions such as Gaussian functions. Instead, it solves the atomic problem numerically and uses the resulting radial functions to construct basis orbitals:
$$
\phi_{i\alpha} \left( \mathbf{r} \right) = R_{\alpha} \left( \left| \mathbf{r} - \mathbf{R}_i \right|  \right) Y_{\ell m} \left( \mathbf{r} - \mathbf{R}_i \right)
$$
centered on atom \(i\). The radial functions \(R_{\alpha}(r)\) are derived from an isolated atom calculation and then truncated at a chosen cutoff radius.

Let's try another equation:
$$F = ma$$

The orbitals used for each atomic species is specified in the OpenMX input file using strings such as `H5.0-s2p1`, where
* a cutoff radius of 5.0 Bohr is specified for hydrogen,
* 2 radial functions are used for the s orbital,
* 1 radial function is used for p orbitals.

Here, we will expect 2 s orbitals and one orbital each for the px, py, and pz orbitals, for a total of 5 orbitals as basis functions. This is essentially a double-zeta s basis plus one set of polarization p functions.

For a system of two H atoms, we might conceptualize the real-space Hamiltonian
elements as $H_{i\alpha,j\beta}$, which is the transition energy between
orbital $\beta$ on the $j$-atom. For Python, we flatten $i\alpha$ to an
integer index using something like this:

```python
# Starting integer index for atom i
i0 = self.orbital_offset[ct_AN]
# Starting integer index for atom j
j0 = self.orbital_offset[Gh_AN]

rows.append(i0 + i) # index orbitals of atom i
cols.append(j0 + j) # index orbitals of atom j
```

Here,
* `ct_AN` maps to atom index $i$. `AN` is for atom number, and `ct` is for central. The central atom is the atom being processed.
* `i` is orbital index $\alpha$ on the central atom,
* `Gh_AN` maps to atom index $j$. `Gh_AN` is for *ghost atom number*. The ghost atom is essentially the neighboring atom whose interaction with the central atom is currently being evaluated.
* `j` is orbital index $\beta$ on the ghost atom.

### Structure of the Basis

OpenMX uses atom-major ordering, then orbital ordering, with index
$\nu = (i,\alpha)$, where

* $i$ indexes atoms, and
* $\alpha$ indexes an orbital on that atom.

Some questions
* What is the atomic problem?

* Break down for me further what it means to use multiple radial functions for
  an s orbital.

* What does it look like to use multiple

## Gauge Convention

## Setup

### Tunneling into Kodiak for Jupyter Service

#### Start a Jupyter Lab server

1. Login on kodiak from a local terminal using something like, but customizing
   `your_username` with a value appropriate to you:

   ```zsh
   ssh your_username@kodiak.baylor.edu
   ```

2. Establish an interactive session on a compute node:

   ```zsh
   qsub -I -q batch
   ```

   This gets us one core on one compute node on the batch queue.

3. Find the hostname for the compute node we're using:

   ```zsh
   hostname
   ```

   For example, I may use `hostname` and see this:

   ```zsh
   $ hostname
   n024
   ```

   Thus, the node is `n024`.

4. Activate the desired Python environment.

   This environment must have `jupyterlab` installed.

5. Start Jupyter (headless)

   ```zsh
   jupyter lab --no-browser --ip=0.0.0.0 --port=8887
   ```

   I'm suggesting that we use a specific (non-default) port number (default
   is 8888) so that we have a distinct Jupyter Lab from the default that may be
   running locally (this likely won't matter).

   The use of `--ip=0.0.0.0` allows the login node to "see" port 8887 of the
   compute node.

   The above command should include in its output information like this:

   > To access the server, open this file in a browser:
   > file:///ion/home/your_username/.local/share/jupyter/runtime/jpserver-1281948-open.html
   >
   > Or copy and paste one of these URLs:
   >
   > http://localhost:8887/lab?token=58ff72ece4f31b1858834096f0dfb92d8897139c0e729438
   >
   > http://127.0.0.1:8887/lab?token=58ff72ece4f31b1858834096f0dfb92d8897139c0e729438

   From this, we save the token. In this example, it's
   `58ff72ece4f31b1858834096f0dfb92d8897139c0e729438`.

6. Tunnel from your laptop or local computer

   ```zsh
   ssh -L 8887:n024:8887 blaire@kodiak.baylor.edu
   ```

   This maps local port 8887 of my laptop to port 8887 of n024 on Kodiak.

7. Point a browser to:

   > http://localhost:8887/lab?token=YOUR_TOKEN

   where `YOUR_TOKEN` is the token we saved earlier.
