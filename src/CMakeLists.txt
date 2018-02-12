# "demo.cxx" and "demo_b.cxx". The extensions are automatically found.
include_directories(inc)

add_executable (dune_sbn dune_sbn.cxx genDune.cxx )

target_link_libraries (dune_sbn SBNfitlib tinyxmllib ${ROOT_LIBRARIES} ${GSL_LIBRARIES})

