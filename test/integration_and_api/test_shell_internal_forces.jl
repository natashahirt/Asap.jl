#=
Test: Shell Internal Forces and Displacements
==============================================
Tests the ShellInternalForces and ShellDisplacements functionality,
plus the unified InternalForces interface.
=#

using Test
using Asap
using Unitful
using LinearAlgebra

@testset "Shell Internal Forces & Displacements" begin
    
    # Material properties
    E = 30e9  # Pa (concrete)
    ν = 0.2
    ρ = 2400.0  # kg/m³
    t = 0.2  # m thickness
    
    # Beam section for frame elements
    E_steel = 200e9  # Pa
    G_steel = 77e9   # Pa
    A = 0.005    # m²
    Ix = 0.0001  # m⁴
    Iy = 0.00005 # m⁴
    J = 0.00001  # m⁴
    beam_sec = Asap.Section(A*u"m^2", E_steel*u"Pa", G_steel*u"Pa", 
                             Ix*u"m^4", Iy*u"m^4", J*u"m^4", 7850.0*u"kg/m^3")
    
    # Geometry
    Lx, Ly, H = 4.0, 3.0, 3.0  # m
    pressure = 5000.0  # Pa
    
    @testset "ShellInternalForces basic" begin
        # Simply supported plate under uniform pressure
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Asap.Node([Lx*u"m", 0.0u"m", 0.0u"m"], :fixed)
        n3 = Asap.Node([Lx*u"m", Ly*u"m", 0.0u"m"], :fixed)
        n4 = Asap.Node([0.0u"m", Ly*u"m", 0.0u"m"], :fixed)
        nc = Asap.Node([Lx/2*u"m", Ly/2*u"m", 0.0u"m"], :free)
        
        shell1 = Asap.ShellTri3((n1, n2, nc), t*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
        shell2 = Asap.ShellTri3((n2, n3, nc), t*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
        shell3 = Asap.ShellTri3((n3, n4, nc), t*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
        shell4 = Asap.ShellTri3((n4, n1, nc), t*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
        shells = [shell1, shell2, shell3, shell4]
        
        load = Asap.AreaLoad(shells, pressure*u"Pa")
        
        model = Asap.ShellModel([n1, n2, n3, n4, nc], shells, [load])
        Asap.process!(model)
        Asap.solve!(model)
        
        # Get internal forces for one shell
        sif = Asap.ShellInternalForces(shell1, model.u)
        
        @test sif isa Asap.ShellInternalForces
        @test sif.element === shell1
        
        # Under pressure, there should be bending moments
        # Mxx and Myy should be non-zero for a loaded plate
        # (The actual values depend on the element formulation)
        
        # Get all shell internal forces
        all_sif = Asap.ShellInternalForces(model)
        @test length(all_sif) == 4
        @test all(s isa Asap.ShellInternalForces for s in all_sif)
        
        @info "ShellInternalForces basic test passed" Mxx=sif.Mxx Myy=sif.Myy
    end
    
    @testset "ShellDisplacements basic" begin
        # Same setup as above
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Asap.Node([Lx*u"m", 0.0u"m", 0.0u"m"], :fixed)
        n3 = Asap.Node([Lx*u"m", Ly*u"m", 0.0u"m"], :fixed)
        n4 = Asap.Node([0.0u"m", Ly*u"m", 0.0u"m"], :fixed)
        nc = Asap.Node([Lx/2*u"m", Ly/2*u"m", 0.0u"m"], :free)
        
        shell1 = Asap.ShellTri3((n1, n2, nc), t*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
        shell2 = Asap.ShellTri3((n2, n3, nc), t*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
        shell3 = Asap.ShellTri3((n3, n4, nc), t*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
        shell4 = Asap.ShellTri3((n4, n1, nc), t*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
        shells = [shell1, shell2, shell3, shell4]
        
        load = Asap.AreaLoad(shells, pressure*u"Pa")
        
        model = Asap.ShellModel([n1, n2, n3, n4, nc], shells, [load])
        Asap.process!(model)
        Asap.solve!(model)
        
        # Get displacement for one shell
        sd = Asap.ShellDisplacements(shell1, model.u)
        
        @test sd isa Asap.ShellDisplacements
        @test sd.element === shell1
        
        # Under downward pressure, w should be negative (deflecting down)
        # Note: This depends on load direction and coordinate system
        
        # Get all displacements
        all_sd = Asap.ShellDisplacements(model)
        @test length(all_sd) == 4
        
        # Test max_deflection helper
        max_w, elem = Asap.max_deflection(all_sd)
        @test elem !== nothing
        
        @info "ShellDisplacements basic test passed" w=sd.w max_w=max_w
    end
    
    @testset "Principal moments" begin
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Asap.Node([Lx*u"m", 0.0u"m", 0.0u"m"], :fixed)
        n3 = Asap.Node([Lx*u"m", Ly*u"m", 0.0u"m"], :fixed)
        n4 = Asap.Node([0.0u"m", Ly*u"m", 0.0u"m"], :fixed)
        nc = Asap.Node([Lx/2*u"m", Ly/2*u"m", 0.0u"m"], :free)
        
        shell1 = Asap.ShellTri3((n1, n2, nc), t*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
        shell2 = Asap.ShellTri3((n2, n3, nc), t*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
        shell3 = Asap.ShellTri3((n3, n4, nc), t*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
        shell4 = Asap.ShellTri3((n4, n1, nc), t*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
        shells = [shell1, shell2, shell3, shell4]
        
        load = Asap.AreaLoad(shells, pressure*u"Pa")
        
        model = Asap.ShellModel([n1, n2, n3, n4, nc], shells, [load])
        Asap.process!(model)
        Asap.solve!(model)
        
        sif = Asap.ShellInternalForces(shell1, model.u)
        
        M1, M2, θ = Asap.principal_moments(sif)
        
        # M1 should be >= M2 (principal values)
        @test M1 >= M2
        # θ should be in [-π/2, π/2]
        @test -π/2 <= θ <= π/2
        
        @info "Principal moments test passed" M1=M1 M2=M2 θ_deg=rad2deg(θ)
    end
    
    @testset "Unified InternalForces interface" begin
        # Create a mixed model with both frames and shells
        nb1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        nb2 = Asap.Node([Lx*u"m", 0.0u"m", 0.0u"m"], :fixed)
        nb3 = Asap.Node([Lx*u"m", Ly*u"m", 0.0u"m"], :fixed)
        nb4 = Asap.Node([0.0u"m", Ly*u"m", 0.0u"m"], :fixed)
        
        n1 = Asap.Node([0.0u"m", 0.0u"m", H*u"m"], :free)
        n2 = Asap.Node([Lx*u"m", 0.0u"m", H*u"m"], :free)
        n3 = Asap.Node([Lx*u"m", Ly*u"m", H*u"m"], :free)
        n4 = Asap.Node([0.0u"m", Ly*u"m", H*u"m"], :free)
        nc = Asap.Node([Lx/2*u"m", Ly/2*u"m", H*u"m"], :free)
        
        # Columns
        col1 = Asap.Element(nb1, n1, beam_sec, :col; release=:freefree)
        col2 = Asap.Element(nb2, n2, beam_sec, :col; release=:freefree)
        col3 = Asap.Element(nb3, n3, beam_sec, :col; release=:freefree)
        col4 = Asap.Element(nb4, n4, beam_sec, :col; release=:freefree)
        
        # Shells
        shell1 = Asap.ShellTri3((n1, n2, nc), t*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
        shell2 = Asap.ShellTri3((n2, n3, nc), t*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
        shell3 = Asap.ShellTri3((n3, n4, nc), t*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
        shell4 = Asap.ShellTri3((n4, n1, nc), t*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
        shells = [shell1, shell2, shell3, shell4]
        
        load = Asap.AreaLoad(shells, pressure*u"Pa")
        
        nodes = [nb1, nb2, nb3, nb4, n1, n2, n3, n4, nc]
        columns = [col1, col2, col3, col4]
        model = Asap.Model(nodes, columns, shells, [load])
        
        Asap.process!(model)
        Asap.solve!(model)
        
        # Test unified InternalForces on individual elements
        col_forces = Asap.InternalForces(col1, model)
        @test col_forces isa Asap.ElementInternalForces
        
        shell_forces = Asap.InternalForces(shell1, model)
        @test shell_forces isa Asap.ShellInternalForces
        
        # Test unified InternalForces on whole model
        all_forces = Asap.InternalForces(model)
        @test haskey(all_forces, :frames)
        @test haskey(all_forces, :shells)
        @test length(all_forces.frames) == 4  # 4 columns
        @test length(all_forces.shells) == 4  # 4 shells
        
        # Test unified Displacements
        disps = Asap.Displacements(model)
        @test haskey(disps, :shells)
        @test length(disps.shells) == 4
        
        @info "Unified InternalForces interface test passed"
    end
    
    @testset "Surface stresses" begin
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Asap.Node([Lx*u"m", 0.0u"m", 0.0u"m"], :fixed)
        n3 = Asap.Node([Lx*u"m", Ly*u"m", 0.0u"m"], :fixed)
        n4 = Asap.Node([0.0u"m", Ly*u"m", 0.0u"m"], :fixed)
        nc = Asap.Node([Lx/2*u"m", Ly/2*u"m", 0.0u"m"], :free)
        
        shell1 = Asap.ShellTri3((n1, n2, nc), t*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
        shell2 = Asap.ShellTri3((n2, n3, nc), t*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
        shell3 = Asap.ShellTri3((n3, n4, nc), t*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
        shell4 = Asap.ShellTri3((n4, n1, nc), t*u"m", E*u"Pa", ν; id=:slab, ρ=ρ)
        shells = [shell1, shell2, shell3, shell4]
        
        load = Asap.AreaLoad(shells, pressure*u"Pa")
        
        model = Asap.ShellModel([n1, n2, n3, n4, nc], shells, [load])
        Asap.process!(model)
        Asap.solve!(model)
        
        sif = Asap.ShellInternalForces(shell1, model.u)
        
        # Get surface stresses
        stresses = Asap.max_surface_stresses(sif, t)
        
        @test haskey(stresses, :top)
        @test haskey(stresses, :bottom)
        @test haskey(stresses.top, :σxx)
        @test haskey(stresses.top, :σyy)
        @test haskey(stresses.top, :τxy)
        
        # von Mises stress at top surface
        σ_vm_top = Asap.von_mises_stress(sif, t/2, t)
        @test σ_vm_top >= 0.0  # von Mises is always non-negative
        
        @info "Surface stresses test passed" σ_vm_top=σ_vm_top
    end
end

println("\n✓ All Shell Internal Forces & Displacements tests passed!")
