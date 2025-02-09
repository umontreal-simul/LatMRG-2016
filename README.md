# LatMRG

## Compiling

### Configuring the Build

LatMRG relies on the
[waf meta build system](https://code.google.com/p/waf/) for configuring and
compiling the software source.
Waf is included in the LatMRG source tree, but it depends on
[Python](http://python.org/download), which must be available on the system
on which LatMRG is to be compiled.

The commands below should work verbatim under Linux and MacOS systems.
**Microsoft Windows** users should replace every instance of `./waf` 
with `C:\Python27\python waf`, assuming that the Python executable
(`python.exe`) was installed under `C:\Python27`, or simply with `python waf`
if the Python installation path is accessible from the system `%PATH%`
environment variable.

Change the current directory to the root directory of the package, for example:

	cd latmrg

if you obtained the source code with the `git` command.
If you obtained the source code from the ZIP archive, the directory should be
named `latmrg-master` instead of `latmrg`.
At the root of the source tree lies the `waf` script, manages the build
process.

LatMRG and LatCommon can be configured to do integer and floating-point
computations with native types of with big integers or real numbers.
This is done by executing:

	./waf configure --ntltypes <NTLTYPES>
	
where `<NTLTYPES>` is one of:

- `LLDD`: to perform integer computations with native `long` integers
  and floating-point computations with native `double` numbers;
- `ZZDD`: to perform integer computations with big integers and
  floating-point computations with native `double` numbers;
- `ZZRR`: to perform integer computations with big integers and
  floating-point computations with arbitrary-precision real numbers.

The following options can also be added to `./waf configure`:

- `--out`: directory in which the files created during the build process will
  be placed;
- `--prefix`: directory under which the LatMRG software will be installed after
  compilation;
- `--boost`: directory under which the Boost libraries are installed;
- `--ntl`: directory under which the NTL library is installed.
- `--testu01`: directory under which TestU01 is installed;
* `--build-docs`: generate the documentation (this requires
  [Doxygen](http://www.stack.nl/~dimitri/doxygen/).

For example, if the header file `vector.h` from the NTL library can be found in
the path `/opt/ntl/include/NTL/vector.h`, you should pass `--ntl
/opt/ntl` to the Waf command.

Try:

	./waf --help

for more options.


It is possible to set the `CXX` environment variable to the path to a specific
C++ compiler to be used to build LatMRG, before running the `waf
configure` command.

The above `waf configure` commands configures `waf` for a minimal build,
without documentation nor code examples.  These can be built by
appending the following options to `waf configure`:

Errors will be reported if required software components cannot be found.  In
that case, you should check the Boost, NTL and TestU01 installation paths.


### Building and Installing

Once everything is configured correctly, the following command will build the
LatMRG library and command-line tool:

	./waf build

If the build process completed without errors, LatMRG can be installed to the
directory specified with the `--prefix` option during the configuration step,
with:

	./waf install


## Running LatMRG

The LatMRG executables can be found in the `bin` subdirectory, under the
installation prefix.
These include:

- `lat<NTLTYPES>`: study lattice properties;
- `seek<NTLTYPES>`: seek parameters for generators;

with `<NTLTYPES>` as specified above.
When `<NTLTYPES>` is `ZZDD`, the additional `maxper`, `findmk` and `findmk2`
programs are found.
Refer to the user guide for further detail.

**NOTE:** Under Windows, the programs have an additional `.exe` extension.


Before executing the LatMRG program, it may be necessary to inform the dynamic
linker where to find the Boost, NTL and TestU01 shared libraries.  Under Linux
this is done by appending the paths to the `LD_LIBRARY_PATH` environment
variable, e.g.,

	export LD_LIBRARY_PATH=/opt/boost/lib:/opt/ntl/lib:/opt/testu01/lib:$LD_LIBRARY_PATH

with a Bash-compatible shell.

**Microsoft Windows** users might need to copy the Boost, NTL and TestU01 DLLs into the
same directory (`$HOME/latmrg/bin`, for example) as the executable programs.

**IMPORTANT:** LatMRG relies on [MIRACL](http://www.miracl.com/) for factorization.  The
environment variable `FACTORHOME` must point to a directory that contains the
MIRACL executable programs.

**IMPORTANT:** Some routines of LatMRG will try to open a file named
`primes.dat`.  It is located in the `data` subdirectory of the source
distribution and must be available in the directory in which the LatMRG
programs are executed (this is inconvenient and could be improved).
