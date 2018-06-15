file(REMOVE_RECURSE
  "ran2_fortran.pdb"
  "ran2_fortran"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/ran2_fortran.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
