# "demo.cxx" and "demo_b.cxx". The extensions are automatically found.

add_executable (dune_sbn dune_sbn.cxx)

target_link_libraries (dune_sbn SBNfitlib tinyxmllib ${ROOT_LIBRARIES})

