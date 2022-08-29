project(Multigrid)
set(ROOT ..)
set(OUTPUT ${ROOT}/star)

add_compile_options(-g -W -Wall -std=c++11)
add_executable(star main.cpp multigrid.cpp)
