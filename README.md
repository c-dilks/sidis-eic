# Largex-eic
General purpose analysis software for SIDIS at the EIC

This repository provides a set of common tools for the analysis of both full and
fast simulations, including the following features:

- kinematics reconstruction methods (e.g., leptonic, hadronic, Jacquet-Blondel,
  etc.)
- calculations of SIDIS variables, such as `PhiH` and `qT`, for single
  particles, as well as jet variables
- application of common set of cuts
- ability to specify arbitrary multi-dimensional binning schemes
- outputs include binned histograms, tables, and other data structures such as
  `TTrees`
- An analysis is primarily driven by macros, used to set up the binning and
  other settings

If you prefer to use your own analysis code, but would still like to make use of
the common tools provided in this repository (e.g., kinematics reconstruction),
this is also possible; you only need to stream the data structure you need, most
likely within the event loop. In this situation, it is recommended you fork the
repository (pull requests are also welcome).


# Setup and Dependencies

## Option 1: Use the Singularity or Docker image

- To minimize setup effort, and provide a consistent development environment, 
  a Singularity image is available, which contains all the dependencies
  pre-built, as well as sample ROOT files
  - First run `container/install.sh` to download and build the Singularity image
    - With no arguments, a usage guide will be printed
    - Default image file location is `container/img/`
    - Note that the image size is about 2 GB
    - Images are hosted on [Docker Hub](https://hub.docker.com/r/cjdilks/largex-eic)
      - (the Docker image is hosted, but Singularity can pull it too)
  - Then run `container/shell.sh` to start a shell in the container
    - This will automatically call `source env.sh` upon shell startup, which
      sets environment variables
  - Proceed with the **Building** section below (just type `make`)

- **Alternatively** if you prefer to use Docker:
  - obtain the image using `docker pull cjdilks/largex-eic:dev`
  - start the container using a standard `docker run` command; you can also use
    the script `container/devscripts/dockerShell.sh`, if you find it useful
    - the Docker image was built assuming a default user ID (UID) of 1000; if your
      UID is different (check with the `id` command), your user name in the container
      may be `I have no name!`, but you should still have read/write permission for
      the present working directory; we have not tested working in this condition,
      due to our preference for Singularity, however suggestions how to improve
      are welcome
    - Docker files are also provided, you can follow `container/devscripts/README.md`
      for instructions how to build your own image (which would allow you to change
      the default UID, or anything else you want)
  - once you are in the Docker container, proceed with the **Building** section below

## Option 2: Setup your own environment

- The other option is to manually set up your environment, by downloading and/or
  building all of the necessary dependencies
- Once you have all the dependencies, proceed with the **Building** section
  below

### Dependencies

- **ROOT**: prefer v6.24.02 or later
- **Delphes**:
  - the analysis is capable of reading `delphes` fast simulation output, and also
    provides a simple wrapper for `delphes` to help keep input `hepmc` and output
    `root` files organized
    - it is not required to use the `delphes` wrapper, but `delphes` libraries are
      needed for the analysis of fast simulation data
  - first, make sure you have a build of `delphes` somewhere, preferably in a
    separate directory
  - set environment variables before doing anything, so this repository knows where your
    `delphes` build is: `source env.sh /path/to/delphes/repository`
    - if you do not specify a path to `delphes` repository, it will use a default
      path given in `env.sh`; it is useful to edit this default path for your own
      convenience
    - it will also symlink `delphes` external code, so analysis macros
      will not complain

## Building

- First make sure environment variables are set by calling `source env.sh`
- Build analysis code with `make`
  - It requires a `root` build as well as `delphes` (see above)
  - All classes are found in the `src/` directory

## Quick Start

- If you're ready to try the software hands-on, follow the [tutorials](tutorial/README.md) in 
  the `tutorial/` directory


# Simulation

## Delphes Fast Simulation

- for convenience, the wrapper script `exeDelphes.sh` is provided, which runs
  `delphes` on a given `hepmc` or `hepmc.gz` file, and sets the output file
  names and the appropriate configuration card
  - configuration cards are stored in the `cards/` directory as a submodule
    - clone this `largex-eic` repository with `--recurse-submodules`, or
      if you already have cloned without submodules, execute
      `git submodule update --init` to obtain them
  - environment must be set first (`source env.sh`)
  - run `exeDelphes.sh` with no arguments for usage guide
  - in the script, you may need to change `exeDelphes` to the proper
    executable, e.g., `DelphesHepMC2` or `DelphesHepMC3`, depending
    on the format of your generator input
  - if reading a gunzipped file (`*.hepmc.gz`), this script will automatically
    stream it through `gunzip`, so there is no need to decompress beforehand
- the output will be a `TTree` stored in a `root` file
  - output files will be placed in `datarec/`
  - input `hepmc(.gz)` files can be kept in `datagen/`

## ATHENA Full Simulation

- TODO 



# Analysis Procedure

After simulation, this repository separates the analysis procedure into two
stages: (1) the *Analysis* stage includes the event loop, which processes either
fast or full simulation output, kinematics reconstruction, and your specified
binning scheme, while (2) the *Post-processing* stage includes histogram
drawing, comparisons, table printouts, and any feature you would like to add

The two stages are driven by macros. Example macros will eventually be added;
for now you can assume any macro named `analysis_*.C` or `postprocess_*.C` are
respective macros for the stages.

- **Note**: most macros stored in this repository must be executed from the
  `largex-eic` top directory, not from within their subdirectory, e.g., run
  `root -b -q tutorial/analysis_template.C`; this is because certain library
  and data directory paths are given as relative paths

## Analysis

### Analysis Macro and Class

- the `Analysis` class is the main class that performs the analysis; it is 
  controlled at the macro level
  - a typical analysis macro must do the following:
    - instantiate `Analysis` (with file names, beam energies, crossing angle)
    - set up bin schemes and bins (arbitrary specification, see below)
    - set any other settings (e.g., a maximum number of events to process,
      useful for quick tests)
    - execute the analysis
    - see `src/Analysis.h` for further documentation in comments
  - the output will be a `root` file, filled with `TObjArray`s of
    histograms
    - each `TObjArray` can be for a different subset of events (bin), e.g.,
      different minimum `y` cuts, so that their histograms can be compared and
      divided; you can open the `root` file in a `TBrowser` to browse the
      histograms
    - the `Histos` class is a container for the histograms, and instances of
      `Histos` will also be streamed to `root` files, along with the binning
      scheme (handled by the `BinSet` class); downstream post processing code
      makes use of these streamed objects, rather than the `TObjArray`s

### Bin Specification

- The bins may be specified arbitrarily, using the `BinSet` and `CutDef` classes
  - see example `analysis_*C` macros
  - `CutDef` can store and apply an arbitrary cut for a single variable, such
    as:
    - ranges: `a<x<b` or `|x-a|<b`
    - minimum or maximum: `x>a` or `x<a`
    - no cut (useful for "full" bins)
  - The set of bins for a variable is defined by `BinSet`, a set of bins
    - These bins can be defined arbitrarily, with the help of the `CutDef`
      class; you can either:
      - Automatically define a set of bins, e.g., `N` bins between `a` and `b`
        - Equal width in linear scale
        - Equal width in log scale (useful for `x` and `Q2`)
        - Any custom `TAxis`
      - Manually define each bin
        - example: specific bins in `z` and `pT`:
          - `|z-0.3|<0.1` and `|pT-0.2|<0.05`
          - `|z-0.7|<0.1` and `|pT-0.5|<0.05`
        - example: 3 different `y` minima:
          - `y>0.05`
          - `y>0.03`
          - `y>0` (no cut)
          - note that the arbitrary specification permits bins to overlap, e.g.,
            an event with `y=0.1` will appear in all three bins
- Multi-dimensional binning
  - Binning in multi-dimensions is allowed, e.g., 3D binning in `x`,`Q2`,`z`
  - See [Adage documentation](doc/adage.md) for more information on how multi-dimensional
    binning is handled, as well as the [Adage syntax reference](doc/syntax.md)
  - Be careful of the curse of dimensionality
    - You can restrict the binning in certain dimensions by taking only diagonal
      elements of a matrix of bins (see `diagonal` settings in `src/Analysis.h`)

### Simple Tree

- The `Analysis` class is capable of producing a simple `TTree`, handled by the
  `SimpleTree` class, which can also be useful for analysis
  - As the name suggests, it is a flat tree with a minimal set of variables,
    specifically needed for asymmetry analysis
  - The tree branches are configured to be compatible with 
    [asymmetry analysis code](https://github.com/c-dilks/dispin),
    built on the [BruFit](https://github.com/dglazier/brufit) framework
  - There is a switch in `Analysis` to enable/disable whether this tree is 
    written


# Post-Processing

### Post-Processing Macro and Class

- results processing is handled by the `PostProcessor` class, which does tasks
  such as printing tables of average values, and drawing ratios of histograms
  - this class is steered by `postprocess_*.C` macros, which includes the
    following:
    - instantiate `Analysis`, needed for bin loops and settings
    - instantiate `PostProcessor`, with the specified `root` file that contains
      output from the analysis macro
    - loop over bins and perform actions
- see `src/PostProcessor.h` and `src/PostProcessor.cxx` for available
  post-processing routines; you are welcome to add your own

### Asymmetry Fitting
- the `SimpleTree` output is compatible with [asymmetry
  code](https://github.com/c-dilks/largex-eic-asym), included here as a
  submodule in `asym/`
  - clone this `largex-eic` repository with `--recurse-submodules`, to get
    `largex-eic-asym` and its main dependency `brufit`
  - follow `asym/README.md`

# Contributions

- This repository is in an early stage of development, so bugs and issues are
  likely
- Contributions are welcome via pull requests and issues reporting; you may also
  find it useful to fork the repository for your own purposes, so that you do
  not have to feel limited by existing code (you can still send pull requests
  from a fork)
- Continuous Integration (CI) will trigger on pull requests, which will build
  and test your contribution; see `Actions` tab for workflows for details
- It is recommended to keep up-to-date with developments by browsing the pull
  requests, issues, and viewing the latest commits by going to the `Insights`
  tab, and clicking `Network` to show the branch topology
