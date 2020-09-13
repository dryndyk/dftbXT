set(sources
  BLASSM/blassm.f
  BLASSM/zblassm.f
  BLASSM/matvec.f
  BLASSM/zmatvec.f
  FORMATS/formats.f
  FORMATS/zformats.f
  FORMATS/unary.f
  FORMATS/zunary.f
  FORMATS/sort.f
  INFO/infofun.f
  INFO/zinfofun.f
  INOUT/inout.f
  INOUT/zinout.f
  ITSOL/zilut.f
  ORDERINGS/ccncommon.f
  ORDERINGS/ccn.f
  ORDERINGS/zccn.f
  ORDERINGS/zcolor.f
  ORDERINGS/zdsepart.f)

add_library(sparskit ${sources})

if(WITH_API OR BUILD_SHARED_LIBS)
  install(TARGETS sparskit DESTINATION "${CMAKE_INSTALL_LIBDIR}" EXPORT dftbplus-targets)
endif()

set(EXPORTED_COMPILED_LIBRARIES ${EXPORTED_COMPILED_LIBRARIES} sparskit PARENT_SCOPE)
