# "demo.cxx" and "demo_b.cxx". The extensions are automatically found.
include_directories(inc)
#include_directories(/usr/local/include/globes)

add_executable (dune_sbn dune_sbn.cxx genDune.cxx )
add_executable (probtest probtest.cxx)
add_executable (biprob biprob.cxx)
add_executable (combined combined.cxx genDune.cxx)

target_link_libraries (dune_sbn SBNfitlib tinyxmllib ${ROOT_LIBRARIES} ${GSL_LIBRARIES} )
target_link_libraries (probtest SBNfitlib tinyxmllib ${ROOT_LIBRARIES} ${GSL_LIBRARIES} )
target_link_libraries (biprob SBNfitlib tinyxmllib ${ROOT_LIBRARIES} ${GSL_LIBRARIES} )
target_link_libraries (combined SBNfitlib tinyxmllib ${ROOT_LIBRARIES} ${GSL_LIBRARIES} )

