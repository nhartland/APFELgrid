autoreconf -i
PREFIX=$(apfel-config --prefix)
./configure --prefi=$PREFIX
make && make check && make install
