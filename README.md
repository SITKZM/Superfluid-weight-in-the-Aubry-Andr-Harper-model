# Dxx_AAH.f90

Superfluid weight (conventional + geometric contributions) of a one-dimensional chain superconductor with a quasiperiodic Aubry–André–Harper (AAH) on-site potential, at zero temperature.

This program **sets the order parameter as a parameter inside the source file**. The gap is therefore an input, not a self-consistent output.

---

## 1. Model

Normal-state Bloch Hamiltonian of the supercell, for a twist (Peierls) phase `k`:

```
H0(k)_{ii}   = -mu + lambda * cos(2*pi*(i-1)/tau_a)          i = 1..N
H0(k)_{i,i+1} = -t * exp(i k),   H0(k)_{i+1,i} = -t * exp(-i k)
```

with periodic boundary conditions (`H0(N,1)`, `H0(1,N)` closed the same way).

The incommensurate golden-ratio modulation is replaced by a rational approximant tau_a:

```
tau = (1 + sqrt(5))/2,   M = nint(N/tau),   tau_a = N/M
```

so that the potential is periodic with period `N` (one supercell) and `k` runs over the reduced Brillouin zone `[-pi/N, pi/N)`. Choosing `N` to be a Fibonacci number makes `tau_a` a good approximant (e.g. `N = 987` → `M = 610`, `tau_a ≈ 1.618033`).

The current operator is `J_x = dH0/dk`:

```
Jx_{i,i+1} = -i t exp(i k),   Jx_{i+1,i} = conjg(-i t exp(i k))
```

`H0(k)` is diagonalized with LAPACK `ZHEEVD`, giving band energies `eps_n(k)` (measured from `mu`, which is already inside `H0`) and eigenvectors; `Jx` is then rotated into that band basis.

## 2. Quantity computed

Uniform s-wave pairing `Delta` is assumed, so the BdG quasiparticle energies in the band basis are

```
E_n = sqrt( eps_n^2 + |Delta|^2 )
```

and the zero-temperature superfluid weight is evaluated as

```
                1      ---   ---            2 |Delta|^2 |<n|J_x|m>|^2
D_xx  =  ------------  >     >     ----------------------------------------
          N * Nk       ---k  --- n,m          E_n E_m ( E_n + E_m )
```

The `k` integral uses the midpoint rule over the reduced BZ: `k = -pi/N + dk*(ik + 1/2)`, `dk = (2*pi/N)/Nk`, `ik = 0 .. Nk-1`.

The sum is split into two parts:

| Quantity   | Definition |
|------------|------------|
| `Dxx_tot`  | full double sum over `n, m` |
| `Dxx_conv` | terms with `eps_n = eps_m` (conventional contribution) |
| `Dxx_geom` | `Dxx_tot - Dxx_conv` (geometric contribution) |

Degeneracy is decided numerically by `|eps_n - eps_m| < tol` with `tol = 1.0d-10`, hardcoded in `calculate_superfluid_weight`.

## 3. Parameters set inside the program

Declared as `parameter` near the top of the source; edit and recompile to change them.

| Name        | Default  | Meaning |
|-------------|----------|---------|
| `delta_mag` | `1.0d-2` | uniform order parameter `|Delta|` |
| `Nk`        | `100`    | number of `k` points (= number of supercells) in the reduced BZ |

## 4. Parameters read from `params.in`

The file is read as the namelist group `input_params` from the current working
directory. Used values:

| Name     | Meaning |
|----------|---------|
| `N`      | number of sites (supercell size); use a Fibonacci number |
| `t`      | nearest-neighbour hopping |
| `mu`     | chemical potential |
| `lambda` | amplitude of the AAH potential |

`temperature`, `U`, `a`, `b`, `n_c`, `criterion`, `Nstop` are **not used** here.
They are kept in the namelist declaration only so that the same `params.in` can be shared with the other programs; a namelist read fails if the file contains a name that is not in the group.

Example `params.in`:

```
&input_params
  temperature = 0.0,
  N = 987,
  U = 0.2d0,
  t = 1.0d0,
  mu = 0.0d0,
  lambda = 0.0d0,
  a = 2.1d0
  b = 0.0d0,
  n_c = 1001,
  criterion = 1e-06,
  Nstop = 1000
/
```

## 5. Build

```
ifx -o b.out Dxx_AAH.f90 -qmkl -qopenmp -mcmodel=medium -shared-intel -heap-arrays
```

`-qopenmp` is only there to enable the threaded MKL; the source itself contains no OpenMP directives.

With gfortran and a reference LAPACK:

```
gfortran -O2 -o b.out Dxx_AAH.f90 -llapack -lblas
```

## 6. Run

```
mkdir -p Dxx_AAH     # the output directory must exist beforehand
./b.out
```

The program prints the progress (`ik+1 completed.`) for each `k` point and the
total run time.

## 7. Output

One file per run:

```
./Dxx_AAH/Delta=<delta_mag>_mu=<mu>_lambda=<lambda>.txt
```

e.g. `./Dxx_AAH/Delta=0.0100_mu=0.00_lambda=0.00.txt`, containing

```
# Dxx_tot        Dxx_conv        Dxx_geom
   1.234567E-03    2.345678E-04    1.000000E-03
```

The name is formatted with `F6.4` for `delta_mag` and `F4.2` for `mu` and `lambda`, so `delta_mag < 10` and `0 <= mu, lambda < 10` are assumed; outside that range the field overflows into `******`. Reruns with the same `delta_mag`, `mu`, `lambda` overwrite the existing file.

## 8. Notes and limitations

- **Zero temperature only.** The `temperature` entry in `params.in` is ignored.
- **`Delta` is not self-consistent.** It is a free input, so a given (`delta_mag`, `mu`, `lambda`) does not correspond to any particular `U` unless you take the value from a self-consistent run (`Self_consistent.f90`).
- **Uniform `Delta` assumed.** Site-dependent pairing, which is what the AAH potential generally produces, is outside the scope of this program — the whole `E_n = sqrt(eps_n^2 + |Delta|^2)` shortcut relies on uniformity.
- **Cost.** `ZHEEVD` and the two `matmul` calls are `O(N^3)` per `k` point, so the run time scales as `O(Nk * N^3)`. Memory is dominated by a few `N x N` complex arrays plus the `ZHEEVD` workspace (`lwork = N^2 + 2N` complex, `lrwork = 2N^2 + 5N + 1` double), i.e. roughly 80 MB at `N = 987`. Build with `-mcmodel=medium -heap-arrays` as in the command above.
- `info` returned by `ZHEEVD` is not checked.
