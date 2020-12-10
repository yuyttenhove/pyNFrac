#! @PYTHON_EXECUTABLE@

################################################################################
# This file is part of pyNFrac
# Copyright (C) 2016 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
#
# pyNFrac is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# pyNFrac is distributed in the hope that it will be useful,
# but WITOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with pyNFrac. If not, see <http://www.gnu.org/licenses/>.
################################################################################

##
# @file test_pyNFrac.py
#
# @brief Unit test for pyNFrac.
#
# @author Bert Vandenbroucke (bv7@st-andrews.ac.uk)
##

import load_module as pyNFrac

import numpy as np
import sys

##
# @brief Unit test for pyNFrac.
##
def main():
  # this file will only be present and contain sensible values if
  # testNeutralFracTable worked
  data = np.loadtxt("testcurve.txt")

  nfracH = pyNFrac.get_neutral_fractions(data[:,0], data[:,1], data[:,2],
                                         data[:,4], data[0,3])

  for i in range(len(nfracH)):
    if abs(nfracH[i] - data[i,5]) > 1.e-5*abs(nfracH[i] + data[i,5]):
      print("Wrong neutral fraction: {val} ({tval}," \
            " relative difference {rdiff})!".format(
        val = nfracH[i], tval = data[i,5],
        rdiff = abs(nfracH[i] - data[i,5]) / abs(nfracH[i] + data[i,5])))
      sys.exit(1)

  sys.exit(0)

##
# @brief Make sure the main function is called.
##
if __name__ == "__main__":
  main()
