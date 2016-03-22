ACLOCAL_AMFLAGS = -I m4

AM_CXXFLAGS = $(LHAPDF_CXXFLAGS) $(ROOT_CXXFLAGS) $(APFEL_CXXFLAGS) $(APPLGRID_CXXFLAGS) $(SIMD_FLAGS)
AM_CPPFLAGS = $(LHAPDF_CPPFLAGS) $(ROOT_CPPFLAGS) $(APFEL_CPPFLAGS) $(APPLGRID_CPPFLAGS) $(SIMD_FLAGS)
AM_LDFLAGS  = $(LHAPDF_LDFLAGS) $(ROOT_LDFLAGS) $(APFEL_LDFLAGS) $(APPLGRID_LDFLAGS) $(SIMD_FLAGS)

lib_LTLIBRARIES = libAPFELgrid.la
libAPFELgrid_la_SOURCES = src/APFELgrid.cc src/fastkernel.cc src/fkgenerator.cc
libAPFELgrid_la_LDFLAGS = -version-info 0:0:0
libAPFELgrid_la_CXXFLAGS = $(AM_CXXFLAGS)
libAPFELgrid_la_CPPFLAGS = $(AM_CPPFLAGS)

check_PROGRAMS = example_gen example_conv
example_gen_SOURCES = tests/example_gen.cc
example_gen_CPPFLAGS = $(AM_CPPFLAGS) -I ./src
example_gen_CXXFLAGS = $(AM_CXXFLAGS) -I ./src
example_gen_LDFLAGS =  $(AM_LDFLAGS) -lAPFELgrid

example_conv_SOURCES = tests/example_conv.cc
example_conv_CPPFLAGS = $(AM_CPPFLAGS) -I ./src
example_conv_CXXFLAGS = $(AM_CXXFLAGS) -I ./src
example_conv_LDFLAGS =  $(AM_LDFLAGS)  -lAPFELgrid

TESTS= tests/fetchTestData.sh $(check_PROGRAMS) tests/clearTestData.sh
