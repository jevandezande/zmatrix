#!/usr/bin/env python3
import sys

import pytest

from pytest import approx

sys.path.insert(0, '..')

from converter import Converter


class TestZmatrix:
    def setup(self):
        self.conv = Converter()
        self.molecules = [
            'methane',
            'ethanol',
            'bromo-chloro-fluoro-methane',
            'Fe2CO9',
            'benzene',
            'cholesterol',
        ]

    def test_read_zmatrix(self):
        conv = self.conv
        lengths = [5, 9, 5, 20, 12, 74]
        for mol, length in zip(self.molecules, lengths):
            conv.read_zmatrix(f'zmatrix/{mol}.dat')
            assert len(conv.zmatrix) == length

    def test_read_cartesian(self):
        conv = self.conv
        lengths = [5, 9, 5, 20, 12, 74]
        for mol, length in zip(self.molecules, lengths):
            conv.read_cartesian(f'cartesian/{mol}.dat')
            assert len(conv.cartesian) == length

    def test_zmatrix_to_cartesian(self):
        conv = self.conv
        for mol in self.molecules:
            conv.read_zmatrix(f'zmatrix/{mol}.dat')
            conv.zmatrix_to_cartesian()
            assert len(conv.zmatrix) == len(conv.cartesian)
            center_of_mass = [0, 0, 0]
            for (z_atom, pos, z_mass), (c_atom, xyz, c_mass) in zip(conv.zmatrix, conv.cartesian):
                assert z_atom == c_atom
                assert z_mass == c_mass
                center_of_mass += c_mass*xyz
            assert center_of_mass == approx([0, 0, 0], abs=1e-10)

    def test_cartesian_to_zmatrix(self):
        conv = self.conv
        for mol in self.molecules:
            conv.read_cartesian(f'cartesian/{mol}.dat')
            conv.cartesian_to_zmatrix()
            assert len(conv.zmatrix) == len(conv.cartesian)
            for i, ((z_atom, pos, z_mass), (c_atom, xyz, c_mass)) in \
                    enumerate(zip(conv.zmatrix, conv.cartesian)):
                assert z_atom == c_atom
                assert z_mass == c_mass
                for connection in pos:
                    if connection:
                        assert connection[0] < i

    def test_str_cartesian(self):
        conv = self.conv
        for mol in self.molecules:
            conv.read_cartesian(f'cartesian/{mol}.dat')
            conv.str_cartesian()
            assert len(string.splitlines()) == len(conv.cartesian)

    def test_str_cartesian(self):
        conv = self.conv
        for mol in self.molecules:
            conv.read_zmatrix(f'zmatrix/{mol}.dat')
            string = conv.str_zmatrix()
            assert len(string.splitlines()) == len(conv.zmatrix)
