#=
Test: AreaLoad API
==================
Tests the unified AreaLoad interface for shell pressure loads.
=#

using Test
using Asap
using Unitful
using LinearAlgebra

@testset "AreaLoad API" begin
    
    # Common setup
    E_concrete = 30e9  # Pa
    ν_concrete = 0.2
    ρ_concrete = 2400.0  # kg/m³
    t_slab = 0.2  # m
    
    E_steel = 200e9  # Pa
    G_steel = 77e9   # Pa
    ρ_steel = 7850.0 # kg/m³
    
    # Beam section (W12x26 equivalent)
    A_beam = 0.005    # m²
    Ix = 0.0001       # m⁴
    Iy = 0.00005      # m⁴
    J = 0.00001       # m⁴
    beam_sec = Asap.Section(A_beam*u"m^2", E_steel*u"Pa", G_steel*u"Pa", 
                             Ix*u"m^4", Iy*u"m^4", J*u"m^4", ρ_steel*u"kg/m^3")
    
    # Column section
    col_sec = beam_sec
    
    Lx, Ly, H = 4.0, 3.0, 3.0  # m
    pressure = 5000.0  # Pa
    total_expected_load = pressure * Lx * Ly  # 60000 N
    
    @testset "AreaLoad with distribute_to=:nodes (FEM)" begin
        # Slab on a platform of columns - each slab node has a column below
        # Column bases are fixed, slab-column connections are pinned (released)
        
        # Column base nodes (fixed)
        nb1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        nb2 = Asap.Node([Lx*u"m", 0.0u"m", 0.0u"m"], :fixed)
        nb3 = Asap.Node([Lx*u"m", Ly*u"m", 0.0u"m"], :fixed)
        nb4 = Asap.Node([0.0u"m", Ly*u"m", 0.0u"m"], :fixed)
        nbc = Asap.Node([Lx/2*u"m", Ly/2*u"m", 0.0u"m"], :fixed)
        
        # Slab nodes (free, supported by columns)
        n1 = Asap.Node([0.0u"m", 0.0u"m", H*u"m"], :free)
        n2 = Asap.Node([Lx*u"m", 0.0u"m", H*u"m"], :free)
        n3 = Asap.Node([Lx*u"m", Ly*u"m", H*u"m"], :free)
        n4 = Asap.Node([0.0u"m", Ly*u"m", H*u"m"], :free)
        nc = Asap.Node([Lx/2*u"m", Ly/2*u"m", H*u"m"], :free)
        
        # Columns (moments released at both ends for pinned behavior)
        col1 = Asap.Element(nb1, n1, col_sec, :col; release=:freefree)
        col2 = Asap.Element(nb2, n2, col_sec, :col; release=:freefree)
        col3 = Asap.Element(nb3, n3, col_sec, :col; release=:freefree)
        col4 = Asap.Element(nb4, n4, col_sec, :col; release=:freefree)
        colc = Asap.Element(nbc, nc, col_sec, :col; release=:freefree)
        
        # Shell panels
        shell1 = Asap.ShellTri3((n1, n2, nc), t_slab*u"m", E_concrete*u"Pa", ν_concrete; id=:slab, ρ=ρ_concrete)
        shell2 = Asap.ShellTri3((n2, n3, nc), t_slab*u"m", E_concrete*u"Pa", ν_concrete; id=:slab, ρ=ρ_concrete)
        shell3 = Asap.ShellTri3((n3, n4, nc), t_slab*u"m", E_concrete*u"Pa", ν_concrete; id=:slab, ρ=ρ_concrete)
        shell4 = Asap.ShellTri3((n4, n1, nc), t_slab*u"m", E_concrete*u"Pa", ν_concrete; id=:slab, ρ=ρ_concrete)
        shells = [shell1, shell2, shell3, shell4]
        
        # Create AreaLoad with default distribute_to=:nodes
        load = Asap.AreaLoad(shells, pressure*u"Pa")
        
        @test load.distribute_to == :nodes
        @test isnothing(load.axis)
        
        # Create mixed model (shells + columns)
        nodes = [nb1, nb2, nb3, nb4, nbc, n1, n2, n3, n4, nc]
        columns = [col1, col2, col3, col4, colc]
        model = Asap.Model(nodes, columns, shells, [load])
        
        Asap.process!(model)
        Asap.solve!(model)
        
        # Check total reaction at column bases equals total load
        R_z = sum(model.reactions[n.globalID[3]] for n in [nb1, nb2, nb3, nb4, nbc])
        @test isapprox(abs(R_z), total_expected_load, rtol=0.01)
        
        @info "AreaLoad FEM test passed" total_reaction=R_z expected=total_expected_load
    end
    
    @testset "AreaLoad with distribute_to=beams (Tributary)" begin
        # Create frame with beams + shells
        n1_base = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2_base = Asap.Node([Lx*u"m", 0.0u"m", 0.0u"m"], :fixed)
        n3_base = Asap.Node([Lx*u"m", Ly*u"m", 0.0u"m"], :fixed)
        n4_base = Asap.Node([0.0u"m", Ly*u"m", 0.0u"m"], :fixed)
        
        n1_top = Asap.Node([0.0u"m", 0.0u"m", H*u"m"], :free)
        n2_top = Asap.Node([Lx*u"m", 0.0u"m", H*u"m"], :free)
        n3_top = Asap.Node([Lx*u"m", Ly*u"m", H*u"m"], :free)
        n4_top = Asap.Node([0.0u"m", Ly*u"m", H*u"m"], :free)
        nc_top = Asap.Node([Lx/2*u"m", Ly/2*u"m", H*u"m"], :free)
        
        # Columns
        col1 = Asap.Element(n1_base, n1_top, col_sec, :col)
        col2 = Asap.Element(n2_base, n2_top, col_sec, :col)
        col3 = Asap.Element(n3_base, n3_top, col_sec, :col)
        col4 = Asap.Element(n4_base, n4_top, col_sec, :col)
        
        # Beams (edge)
        beam1 = Asap.Element(n1_top, n2_top, beam_sec, :beam)
        beam2 = Asap.Element(n2_top, n3_top, beam_sec, :beam)
        beam3 = Asap.Element(n3_top, n4_top, beam_sec, :beam)
        beam4 = Asap.Element(n4_top, n1_top, beam_sec, :beam)
        beams = [beam1, beam2, beam3, beam4]
        
        # Shells (for geometry only in tributary mode)
        shell1 = Asap.ShellTri3((n1_top, n2_top, nc_top), t_slab*u"m", E_concrete*u"Pa", ν_concrete; id=:slab, ρ=ρ_concrete)
        shell2 = Asap.ShellTri3((n2_top, n3_top, nc_top), t_slab*u"m", E_concrete*u"Pa", ν_concrete; id=:slab, ρ=ρ_concrete)
        shell3 = Asap.ShellTri3((n3_top, n4_top, nc_top), t_slab*u"m", E_concrete*u"Pa", ν_concrete; id=:slab, ρ=ρ_concrete)
        shell4 = Asap.ShellTri3((n4_top, n1_top, nc_top), t_slab*u"m", E_concrete*u"Pa", ν_concrete; id=:slab, ρ=ρ_concrete)
        shells = [shell1, shell2, shell3, shell4]
        
        # Create AreaLoad with distribute_to=beams
        load = Asap.AreaLoad(shells, pressure*u"Pa"; distribute_to=beams)
        
        @test load.distribute_to === beams
        @test isnothing(load.axis)
        
        # Create frame model (no shells in model, shells only define geometry for tributary)
        nodes = [n1_base, n2_base, n3_base, n4_base, n1_top, n2_top, n3_top, n4_top]
        elements = [col1, col2, col3, col4, beam1, beam2, beam3, beam4]
        model = Asap.FrameModel(nodes, elements, [load])
        
        Asap.process!(model)
        Asap.solve!(model)
        Asap.post_process!(model)
        
        # Check total reaction
        R_z = sum(model.reactions[n.globalID[3]] for n in [n1_base, n2_base, n3_base, n4_base])
        @test isapprox(abs(R_z), total_expected_load, rtol=0.01)
        
        # Check beams have shear (load went through beams, not bypassed)
        Vy_beam1 = abs(beam1.forces[2])  # Local shear
        @test Vy_beam1 > 1000.0  # Significant shear
        
        @info "AreaLoad Tributary test passed" total_reaction=R_z beam_shear=Vy_beam1
    end
    
    @testset "AreaLoad with one-way axis" begin
        # Same setup as above
        n1_base = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2_base = Asap.Node([Lx*u"m", 0.0u"m", 0.0u"m"], :fixed)
        n3_base = Asap.Node([Lx*u"m", Ly*u"m", 0.0u"m"], :fixed)
        n4_base = Asap.Node([0.0u"m", Ly*u"m", 0.0u"m"], :fixed)
        
        n1_top = Asap.Node([0.0u"m", 0.0u"m", H*u"m"], :free)
        n2_top = Asap.Node([Lx*u"m", 0.0u"m", H*u"m"], :free)
        n3_top = Asap.Node([Lx*u"m", Ly*u"m", H*u"m"], :free)
        n4_top = Asap.Node([0.0u"m", Ly*u"m", H*u"m"], :free)
        nc_top = Asap.Node([Lx/2*u"m", Ly/2*u"m", H*u"m"], :free)
        
        col1 = Asap.Element(n1_base, n1_top, col_sec, :col)
        col2 = Asap.Element(n2_base, n2_top, col_sec, :col)
        col3 = Asap.Element(n3_base, n3_top, col_sec, :col)
        col4 = Asap.Element(n4_base, n4_top, col_sec, :col)
        
        beam1 = Asap.Element(n1_top, n2_top, beam_sec, :beam_x1)  # X-parallel (no load for axis=(1,0))
        beam2 = Asap.Element(n2_top, n3_top, beam_sec, :beam_y1)  # Y-parallel (gets load)
        beam3 = Asap.Element(n3_top, n4_top, beam_sec, :beam_x2)  # X-parallel (no load)
        beam4 = Asap.Element(n4_top, n1_top, beam_sec, :beam_y2)  # Y-parallel (gets load)
        beams = [beam1, beam2, beam3, beam4]
        
        shell1 = Asap.ShellTri3((n1_top, n2_top, nc_top), t_slab*u"m", E_concrete*u"Pa", ν_concrete; id=:slab, ρ=ρ_concrete)
        shell2 = Asap.ShellTri3((n2_top, n3_top, nc_top), t_slab*u"m", E_concrete*u"Pa", ν_concrete; id=:slab, ρ=ρ_concrete)
        shell3 = Asap.ShellTri3((n3_top, n4_top, nc_top), t_slab*u"m", E_concrete*u"Pa", ν_concrete; id=:slab, ρ=ρ_concrete)
        shell4 = Asap.ShellTri3((n4_top, n1_top, nc_top), t_slab*u"m", E_concrete*u"Pa", ν_concrete; id=:slab, ρ=ρ_concrete)
        shells = [shell1, shell2, shell3, shell4]
        
        # One-way along X: load goes to Y-parallel edges (beams 2 and 4)
        load = Asap.AreaLoad(shells, pressure*u"Pa"; distribute_to=beams, axis=(1.0, 0.0))
        
        @test load.axis == (1.0, 0.0)
        
        nodes = [n1_base, n2_base, n3_base, n4_base, n1_top, n2_top, n3_top, n4_top]
        elements = [col1, col2, col3, col4, beam1, beam2, beam3, beam4]
        model = Asap.FrameModel(nodes, elements, [load])
        
        Asap.process!(model)
        Asap.solve!(model)
        Asap.post_process!(model)
        
        # Total reaction should still equal total load
        R_z = sum(model.reactions[n.globalID[3]] for n in [n1_base, n2_base, n3_base, n4_base])
        @test isapprox(abs(R_z), total_expected_load, rtol=0.01)
        
        # Y-parallel beams should have significant shear
        Vy_beam2 = abs(beam2.forces[2])
        Vy_beam4 = abs(beam4.forces[2])
        
        # X-parallel beams should have minimal shear (near zero)
        Vy_beam1 = abs(beam1.forces[2])
        Vy_beam3 = abs(beam3.forces[2])
        
        @test Vy_beam2 > 5000.0  # Y-beam gets load
        @test Vy_beam4 > 5000.0  # Y-beam gets load
        @test Vy_beam1 < 100.0   # X-beam minimal
        @test Vy_beam3 < 100.0   # X-beam minimal
        
        @info "AreaLoad One-Way test passed" Vy_Y_beams=(Vy_beam2, Vy_beam4) Vy_X_beams=(Vy_beam1, Vy_beam3)
    end
    
    @testset "AreaLoad single shell (backward compat)" begin
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Asap.Node([Lx*u"m", 0.0u"m", 0.0u"m"], :fixed)
        n3 = Asap.Node([Lx*u"m", Ly*u"m", 0.0u"m"], :fixed)
        
        shell = Asap.ShellTri3((n1, n2, n3), t_slab*u"m", E_concrete*u"Pa", ν_concrete; id=:slab, ρ=ρ_concrete)
        
        # Single shell (not vector)
        load = Asap.AreaLoad(shell, pressure*u"Pa")
        
        @test length(load.shells) == 1
        @test load.shells[1] === shell
        
        @info "AreaLoad single-shell test passed"
    end
    
    @testset "SurfaceLoad backward compatibility" begin
        n1 = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2 = Asap.Node([Lx*u"m", 0.0u"m", 0.0u"m"], :fixed)
        n3 = Asap.Node([Lx*u"m", Ly*u"m", 0.0u"m"], :fixed)
        
        shell = Asap.ShellTri3((n1, n2, n3), t_slab*u"m", E_concrete*u"Pa", ν_concrete; id=:slab, ρ=ρ_concrete)
        
        # Old SurfaceLoad API should still work
        load = Asap.SurfaceLoad(shell, pressure*u"Pa")
        
        @test load isa Asap.AreaLoad
        @test load.distribute_to == :nodes
        
        @info "SurfaceLoad backward compat test passed"
    end
    
    @testset "Mixed Model with AreaLoad (FEM mode)" begin
        # Full mixed model: shells for stiffness + beams
        n1_base = Asap.Node([0.0u"m", 0.0u"m", 0.0u"m"], :fixed)
        n2_base = Asap.Node([Lx*u"m", 0.0u"m", 0.0u"m"], :fixed)
        n3_base = Asap.Node([Lx*u"m", Ly*u"m", 0.0u"m"], :fixed)
        n4_base = Asap.Node([0.0u"m", Ly*u"m", 0.0u"m"], :fixed)
        
        n1_top = Asap.Node([0.0u"m", 0.0u"m", H*u"m"], :free)
        n2_top = Asap.Node([Lx*u"m", 0.0u"m", H*u"m"], :free)
        n3_top = Asap.Node([Lx*u"m", Ly*u"m", H*u"m"], :free)
        n4_top = Asap.Node([0.0u"m", Ly*u"m", H*u"m"], :free)
        nc_top = Asap.Node([Lx/2*u"m", Ly/2*u"m", H*u"m"], :free)
        
        col1 = Asap.Element(n1_base, n1_top, col_sec, :col)
        col2 = Asap.Element(n2_base, n2_top, col_sec, :col)
        col3 = Asap.Element(n3_base, n3_top, col_sec, :col)
        col4 = Asap.Element(n4_base, n4_top, col_sec, :col)
        
        beam1 = Asap.Element(n1_top, n2_top, beam_sec, :beam)
        beam2 = Asap.Element(n2_top, n3_top, beam_sec, :beam)
        beam3 = Asap.Element(n3_top, n4_top, beam_sec, :beam)
        beam4 = Asap.Element(n4_top, n1_top, beam_sec, :beam)
        
        shell1 = Asap.ShellTri3((n1_top, n2_top, nc_top), t_slab*u"m", E_concrete*u"Pa", ν_concrete; id=:slab, ρ=ρ_concrete)
        shell2 = Asap.ShellTri3((n2_top, n3_top, nc_top), t_slab*u"m", E_concrete*u"Pa", ν_concrete; id=:slab, ρ=ρ_concrete)
        shell3 = Asap.ShellTri3((n3_top, n4_top, nc_top), t_slab*u"m", E_concrete*u"Pa", ν_concrete; id=:slab, ρ=ρ_concrete)
        shell4 = Asap.ShellTri3((n4_top, n1_top, nc_top), t_slab*u"m", E_concrete*u"Pa", ν_concrete; id=:slab, ρ=ρ_concrete)
        shells = [shell1, shell2, shell3, shell4]
        
        # AreaLoad in FEM mode on mixed model
        load = Asap.AreaLoad(shells, pressure*u"Pa")  # distribute_to=:nodes
        
        nodes = [n1_base, n2_base, n3_base, n4_base, n1_top, n2_top, n3_top, n4_top, nc_top]
        frame_elements = [col1, col2, col3, col4, beam1, beam2, beam3, beam4]
        
        # Mixed model with both shells and frames
        model = Asap.Model(nodes, frame_elements, shells, [load])
        
        @test Asap.has_frame_elements(model)
        @test Asap.has_shell_elements(model)
        @test Asap.is_mixed(model)
        
        Asap.process!(model)
        Asap.solve!(model)
        
        # Check equilibrium
        R_z = sum(model.reactions[n.globalID[3]] for n in [n1_base, n2_base, n3_base, n4_base])
        @test isapprox(abs(R_z), total_expected_load, rtol=0.01)
        
        @info "Mixed Model with AreaLoad test passed" total_reaction=R_z
    end
end

println("\n✓ All AreaLoad API tests passed!")
