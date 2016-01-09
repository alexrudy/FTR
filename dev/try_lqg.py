#!/usr/bin/env python
# -*- coding: utf-8 -*-

from FTR import lqg

f = lqg.LQGFilter.generate_integrator(0.5, 0.995, (10,10))
f.write("LQG-Test.hdf5")

f2 = lqg.LQGFilter.read("LQG-Test.hdf5")
print(f2)

