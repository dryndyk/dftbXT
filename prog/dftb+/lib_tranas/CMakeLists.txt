set(curdir "lib_tranas")

set(sources-fpp
  ${curdir}/clock.F90
  ${curdir}/complexbands.F90
  ${curdir}/contselfenergy.F90
  ${curdir}/distributions.F90
  ${curdir}/elphdb.F90
  ${curdir}/elphdd.F90
  ${curdir}/elphds.F90
  ${curdir}/elph.F90
  ${curdir}/energy_mesh.F90
  ${curdir}/globals.F90
  ${curdir}/input_output.F90
  ${curdir}/inversions.F90
  ${curdir}/lib_param.F90
  ${curdir}/ln_allocation.F90
  ${curdir}/ln_constants.F90
  ${curdir}/ln_extract.F90
  ${curdir}/ln_precision.F90
  ${curdir}/ln_structure.F90
  ${curdir}/load.F90
  ${curdir}/lowdin.F90
  ${curdir}/mat_def.F90
  ${curdir}/metis_interface.F90
  ${curdir}/mkl_dfti.F90
  ${curdir}/mpi_globals.F90
  ${curdir}/outmatrix.F90
  ${curdir}/phph.F90
  ${curdir}/population.F90
  ${curdir}/rcm_module.F90
  ${curdir}/scattstates.F90
  ${curdir}/sparsekit_drv.F90
  ${curdir}/tranas.F90
  ${curdir}/tranas_ngf_integrations.F90
  ${curdir}/tranas_ngf_iterative.F90
  ${curdir}/tranas_ngf_mbngf.F90
  ${curdir}/tranas_types_main.F90
  ${curdir}/tranas_types_mbngf.F90
  ${curdir}/transform.F90)

set(ALL-SOURCES-FPP ${ALL-SOURCES-FPP} ${sources-fpp} PARENT_SCOPE)

