#!/usr/bin/python

from converter import Converter

a = Converter()
a.run_zmatrix( 'tests/Fe2CO9-zmatrix.dat', 'tests/Fe2CO9-cartesian-out.dat' )
#a.run_cartesian( 'tests/Fe2CO9-cartesian.dat', 'tests/Fe2CO9-zmatrix-out.dat' )
