---
title: "Continuous integration"
authors:
- Guido Cossu
- Ed Bennett
---

Every push of the develop branch of Grid,
and every pull request,
is tested using [TeamCity][teamcity].
The [deployed TeamCity instance][grid-teamcity]
is attached to the [Tursa][tursa] facility,
where workflow runs take place;
guest login is available to view build statistics.

The current status is: 
[![Teamcity status](https://ci.dev.dirac.ed.ac.uk/guestAuth/app/rest/builds/aggregated/strob:(buildType:(affectedProject(id:GridBasedSoftware_Grid)),branch:default:true)/statusIcon.svg)](https://ci.dev.dirac.ed.ac.uk/project/GridBasedSoftware_Grid?mode=builds) 

## CI tasks

### Build and test

TeamCity will test the compilation workflow for SU(3),
for single and double precision versions,
on the following compilers:

- Clang 10 targeting AVX2
- NVCC 11 targeting NVIDIA A100

The system runs a subset of the test suite on CPU.

### Benchmark

TeamCity also runs benchmarks on both CPU and GPU,
with the following node types:

- 2 &times; AMD EPYC 48-core
- 4 &times; NVIDIA A100-40GB
- 4 &times; NVIDIA A100-80GB

In each case,
tests run both on a single node,
and on two nodes with MPI,
verifying both single-node and multi-node performance.
Benchmark statistics are retained,
allowing to view performance improvements/regressions over time,
both in the software and the machine configuration.

### Continuous deployment

Builds from TeamCity are available to compile against for Tursa users.
These are provided as a [Pixi][pixi] environment,
under `/mnt/lustre/tursafs1/home/y07/shared/grid/`.

[grid-teamcity]: <https://ci.dev.dirac.ed.ac.uk/>
[pixi]: <https://pixi.prefix.dev/latest/>
[teamcity]: <https://www.jetbrains.com/teamcity/specials/teamcity/teamcity.html>
[tursa]: <https://www.epcc.ed.ac.uk/hpc-services/dirac-tursa-gpu>
