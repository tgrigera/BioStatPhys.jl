# test_rdf.jl -- test radial distribution function
#
# Copyright (C) 2024 Tomas S. Grigera <tgrigera@iflysib.unlp.edu.ar>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For details, see the file LICENSE in the root directory, or
# check <https://www.gnu.org/licenses/>.

using BioStatPhys
using StaticArrays
import Statistics
using Random

function test_rdf_pbc_2d(nconf,npart)
    L = 10.
    region = PeriodicRectangle(L,L)
    pos = [MVector(L*rand(),L*rand()) for _=1:npart]
    dc = density_correlation(region,pos,Δr=0.1)
    for _=1:nconf-1
        for i=1:npart pos[i] .= [L*rand(),L*rand()] end
        density_correlation!(dc,pos)
    end
    return rdf(dc)
end

function test_rdf_pbc_2d_rho_vs_gr(rho=true)
    L = 10.
    region = PeriodicRectangle(L,L)
    dc = density_correlation(region,Δr=0.1)
    rcm,_ = rdf(dc)
    rcm = zeros(length((rcm)))
    nparts = [500,700,1000,1500]
    Random.seed!(9904)
    for npart in nparts
        pos = [MVector(L*rand(),L*rand()) for _=1:npart]
        density_correlation!(dc,pos)
        dcs = density_correlation(region,pos,Δr=0.1)
        rr,_ = rdf(dcs,two_particle_density=rho)
        for (i,r) in enumerate(rr) rcm[i]+=r end
    end
    rc,_ = rdf(dc,two_particle_density=rho)
    return rc, rcm/length(nparts)
end

function test_rdf_npbc_2d(nconf,npart)
    L = 10.
    region = Rectangle(L,L;x0=SVector(0.,0.))
    pos = [MVector(L*rand(),L*rand()) for _=1:npart]
    dc = density_correlation(region,pos,Δr=0.1)
    for _=1:nconf-1
        for i=1:npart pos[i] .= [L*rand(),L*rand()] end
        density_correlation!(dc,pos)
    end
    return rdf(dc)
end

function test_rdf_pbc_2d(nconf,npart_min,npart_max)
    L = 10.
    region = PeriodicRectangle(L,L)
    np = rand(npart_min:npart_max)
    pos = [MVector(L*rand(),L*rand()) for _=1:np]
    dc = BioStatPhys.density_correlation(region,pos,Δr=0.1)
    for _=1:nconf-1
        np = rand(npart_min:npart_max)
        pos = [MVector(L*rand(),L*rand()) for _=1:np]
        BioStatPhys.density_correlation!(dc,pos)
    end
    return BioStatPhys.rdf(dc)
end

function test_rdf_npbc_2d(nconf,npart_min,npart_max)
    L = 10.
    region = Rectangle(L,L,x0=SVector(0.,0.))
    np = rand(npart_min:npart_max)
    pos = [MVector(L*rand(),L*rand()) for _=1:np]
    dc = BioStatPhys.density_correlation(region,pos,Δr=0.1)
    for _=1:nconf-1
        np = rand(npart_min:npart_max)
        pos = [MVector(L*rand(),L*rand()) for _=1:np]
        BioStatPhys.density_correlation!(dc,pos)
    end
    return BioStatPhys.rdf(dc)
end
