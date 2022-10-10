# binvectest.jl
#
# Copyright (C) 2022 Tomas S. Grigera <tgrigera@iflysib.unlp.edu.ar>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# For details see the file LICENSE in the root directory, or check
# <https://www.gnu.org/licenses/>.

BinnedVector_expected=Int[11,5,2,2,2,2,2,2,2,2,2,2]

function BinnedVector_test()
    A = BinnedVector{Int}(10,min=0.,max=10.,init=zeros)
    for x=0:0.5:10. A[x]+=1  end
    A[-100.]+=5
    A[100.]+=10
    return A.data==BinnedVector_expected
end
