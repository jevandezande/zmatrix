#!/usr/bin/python

from converter import Converter

a = Converter()
a.run_zmatrix('tests/zmatrix/methane.dat',
              'tests/output/methane-cartesian.dat')
a.run_zmatrix('tests/zmatrix/ethanol.dat',
              'tests/output/ethanol-cartesian.dat')
a.run_zmatrix('tests/zmatrix/bromo-chloro-fluoro-methane.dat',
              'tests/output/bromo-chloro-fluoro-methane-cartesian.dat')
a.run_zmatrix('tests/zmatrix/Fe2CO9.dat', 'tests/output/Fe2CO9-cartesian.dat')
a.run_zmatrix('tests/zmatrix/benzene.dat',
              'tests/output/benzene-cartesian.dat')
a.run_zmatrix('tests/zmatrix/cholesterol.dat',
              'tests/output/cholesterol-cartesian.dat')

a.run_cartesian('tests/cartesian/methane.dat',
                'tests/output/methane-zmatrix.dat')
a.run_cartesian('tests/cartesian/bromo-chloro-fluoro-methane.dat',
                'tests/output/bromo-chloro-fluoro-methane-zmatrix.dat')
a.run_cartesian('tests/cartesian/ethanol.dat',
                'tests/output/ethanol-zmatrix.dat')
a.run_cartesian('tests/cartesian/Fe2CO9.dat',
                'tests/output/Fe2CO9-zmatrix.dat')
a.run_cartesian('tests/cartesian/benzene.dat',
                'tests/output/benzene-zmatrix.dat')
a.run_cartesian('tests/cartesian/cholesterol.dat',
                'tests/output/cholesterol-zmatrix.dat')
