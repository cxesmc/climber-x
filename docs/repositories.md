
# CLIMBER-X repository overview

Repositories are organized into several communities.

- [github.com/cxesmc](https://github.com/cxesmc): the main CLIMBER-X ESM community - for CLIMBER-X specific development.
- [gitlab.pik-potsdam.de/cxesmc](https://gitlab.pik-potsdam.de/cxesmc): necessary for housing large input files (otherwise GitHub community would be enough).
- [github.com/fesmc](https://github.com/fesmc): the Fast ESM community - for shared resources.
- [github.com/palma-ice](https://github.com/palma-ice): the PalMA ice-sheet modeling group at the Complutense University of Madrid, where the ice-sheet model Yelmo is developed.

## Specific repository list

| Repository | Level   | Location    |
| --- | --- | --- |
| Main CLIMBER-X code base | base | [https://github.com/cxesmc/climber-x](https://github.com/cxesmc/climber-x) |
| CLIMBER-X input files | base | [https://gitlab.pik-potsdam.de/cxesmc/climber-x-input](https://gitlab.pik-potsdam.de/cxesmc/climber-x-input) |
| FESM Utilities | base | [https://github.com/fesmc/fesm-utils](https://github.com/fesmc/fesm-utils) |
| Coordinates package¹ | base | [https://github.com/fesmc/coordinates](https://github.com/fesmc/coordinates) |
| Biogeochemistry (BGC, i.e. HAMOCC)² | optional | [https://github.com/cxesmc/bgc](https://github.com/cxesmc/bgc) |
| Solid Earth model (VILMA)³ | optional | [https://github.com/cxesmc/vilma](https://github.com/cxesmc/vilma) |
| Ice-sheet model (Yelmo) | optional | [https://github.com/palma-ice/yelmo](https://github.com/palma-ice/yelmo) |

¹Repository will eventually be merged into `fesm-utils`.

²PRIVATE. Since the HAMOCC model code is not open source, the `bgc` repository is private at the moment and
you need to be given permission in order to access it. HAMOCC is covered by the Max Planck Institute for
Meteorology software licence agreement as part of the MPI-ESM ([https://code.mpimet.mpg.de/attachments/download/26986/MPI-ESM_SLA_v3.4.pdf](https://code.mpimet.mpg.de/attachments/download/26986/MPI-ESM_SLA_v3.4.pdf)).
A pre-requisite to access the `bgc` repository is therefore that you agree to the MPI-ESM license
by following the steps outlined here: [https://code.mpimet.mpg.de/projects/mpi-esm-license](https://code.mpimet.mpg.de/projects/mpi-esm-license).
Once you have done so, send an email to [Matteo Willeit](mailto:matteo.willeit@gmail.com?subject=[GitHub]%20bgc%20source%20code) and you will be granted permission to access the `bgc` repository.

³PRIVATE. Since the VILMA model code is not open source, the `vilma` repository is private at the moment
and you need to be given permission in order to access it. Please send an email
to [Matteo Willeit and Volker Klemann](mailto:matteo.willeit@gmail.com,volker.klemann@gfz.de?subject=[GitHub]%20VILMA%20access)
and you will be granted permission to access the `vilma` repository.
