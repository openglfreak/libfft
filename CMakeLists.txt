cmake_minimum_required(VERSION 2.6)
project(tests)

include_directories(${PROJECT_SOURCE_DIR}/include)

if(0)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -pg -Wall -fno-rtti -std=gnu++11") #
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg -g -Og -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -m64 -flto -fexpensive-optimizations -finline-functions -fno-exceptions -pipe -static -static-libgcc -static-libstdc++") # -fomit-frame-pointer -s
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -g -O3 -fexpensive-optimizations")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -flto")
else()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fno-rtti -std=gnu++11")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -pg -g -Og -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -m64 -flto -fomit-frame-pointer -fexpensive-optimizations -finline-functions -fno-exceptions -pipe -static -static-libgcc -static-libstdc++ -s")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -g -O3 -fexpensive-optimizations")
set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} -flto")
endif()
add_executable(tests ${PROJECT_SOURCE_DIR}/src/tests/main.cpp)
