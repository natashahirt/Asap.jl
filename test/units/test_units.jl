# =============================================================================
# Unit System Tests
# =============================================================================
# Tests for Asap's unit definitions, type aliases, and conversion helpers.
# Asap is the canonical source for structural engineering units in the ecosystem.

@testset "Unit Exports" begin
    # US Customary units should be exported and usable
    @testset "US Customary Units" begin
        # kip (kilopound-force)
        @test 1.0kip ≈ 4448.2216152605u"N" rtol=1e-6
        @test uconvert(u"lbf", 1.0kip) ≈ 1000.0u"lbf" rtol=1e-6
        
        # ksi (kips per square inch)
        @test 1.0ksi ≈ 6.894757e6u"Pa" rtol=1e-6
        @test uconvert(u"psi", 1.0ksi) ≈ 1000.0u"psi" rtol=1e-6
        
        # psf (pounds per square foot)
        @test 1.0psf ≈ 47.88025898u"Pa" rtol=1e-6
        
        # ksf (kips per square foot)
        @test uconvert(psf, 1.0ksf) ≈ 1000.0psf rtol=1e-6
        @test 1.0ksf ≈ 47880.25898u"Pa" rtol=1e-6
        
        # pcf (pounds per cubic foot)
        @test 1.0pcf ≈ 16.01846337u"kg/m^3" rtol=1e-6
    end
    
    # Physical constants
    @testset "Physical Constants" begin
        @test GRAVITY ≈ 9.80665u"m/s^2"
    end
end

@testset "Type Aliases" begin
    # Dimension-based type aliases (for function signatures)
    @testset "Dimension-Based Aliases" begin
        # These should accept any compatible unit
        @test 1.0u"m" isa Length
        @test 1.0u"ft" isa Length
        @test 1.0u"inch" isa Length
        
        @test 1.0u"m^2" isa Area
        @test 1.0u"ft^2" isa Area
        
        @test 1.0u"m^3" isa Volume
        @test 1.0u"ft^3" isa Volume
        @test 1.0u"m^3" isa SectionModulus  # Alias for Volume (L³)
        
        @test 1.0u"Pa" isa Pressure
        @test 1.0ksi isa Pressure
        @test 1.0psf isa Pressure
        
        @test 1.0u"N" isa Force
        @test 1.0kip isa Force
        @test 1.0u"lbf" isa Force
        
        @test 1.0u"N*m" isa Moment
        @test 1.0kip*u"ft" isa Moment
        
        @test 1.0u"N/m" isa LinearLoad
        @test 1.0kip/u"ft" isa LinearLoad
        
        @test 1.0u"kg/m^3" isa Density
        @test 1.0pcf isa Density
        
        @test 1.0u"m/s^2" isa Acceleration
    end
    
    # Concrete type aliases (for struct fields - Float64 precision)
    @testset "Concrete Type Aliases" begin
        @test 1.0u"m" isa LengthQuantity
        @test 1.0u"m^2" isa AreaQuantity
        @test 1.0u"m^3" isa VolumeQuantity
        @test 1.0u"Pa" isa PressureQuantity
        @test 1.0u"N" isa ForceQuantity
        @test 1.0u"N*m" isa MomentQuantity
        @test 1.0u"N/m" isa ForcePerLength
    end
    
    # Section property type aliases
    @testset "Section Property Aliases" begin
        @test 1.0u"m^4" isa SecondMomentOfArea
        @test 1.0u"m^4" isa TorsionalConstant
        @test 1.0u"m^4" isa MomentOfInertia  # Alias
        @test 1.0u"m^6" isa WarpingConstant
    end
end

@testset "Conversion Helpers - US Customary" begin
    # to_inches
    @test to_inches(1.0u"m") ≈ 39.3701 rtol=1e-4
    @test to_inches(12.0u"ft") ≈ 144.0 rtol=1e-6
    @test to_inches(1.0u"inch") ≈ 1.0
    @test to_inches(100.0) == 100.0  # Pass-through for Real
    
    # to_sqinches
    @test to_sqinches(1.0u"m^2") ≈ 1550.0031 rtol=1e-4
    @test to_sqinches(1.0u"ft^2") ≈ 144.0 rtol=1e-6
    @test to_sqinches(100.0) == 100.0  # Pass-through for Real
    
    # to_ksi
    @test to_ksi(1.0ksi) ≈ 1.0
    @test to_ksi(1000.0u"psi") ≈ 1.0 rtol=1e-6
    @test to_ksi(6.894757u"MPa") ≈ 1.0 rtol=1e-5
    @test to_ksi(50.0) == 50.0  # Pass-through for Real
    
    # to_kip
    @test to_kip(1.0kip) ≈ 1.0
    @test to_kip(1000.0u"lbf") ≈ 1.0
    @test to_kip(4.448222u"kN") ≈ 1.0 rtol=1e-4
    @test to_kip(100.0) == 100.0  # Pass-through for Real
    
    # to_kipft
    @test to_kipft(1.0kip*u"ft") ≈ 1.0
    @test to_kipft(12.0kip*u"inch") ≈ 1.0
    @test to_kipft(100.0) == 100.0  # Pass-through for Real
end

@testset "Conversion Helpers - SI" begin
    # to_meters
    @test to_meters(1.0u"m") ≈ 1.0
    @test to_meters(1.0u"ft") ≈ 0.3048 rtol=1e-6
    @test to_meters(1000.0u"mm") ≈ 1.0
    @test to_meters(100.0) == 100.0  # Pass-through for Real
    
    # to_meters_squared
    @test to_meters_squared(1.0u"m^2") ≈ 1.0
    @test to_meters_squared(1.0u"ft^2") ≈ 0.092903 rtol=1e-4
    
    # to_meters_fourth
    @test to_meters_fourth(1.0u"m^4") ≈ 1.0
    @test to_meters_fourth(1.0u"cm^4") ≈ 1e-8 rtol=1e-10
    
    # to_pascals
    @test to_pascals(1.0u"Pa") ≈ 1.0
    @test to_pascals(1.0u"kPa") ≈ 1000.0
    @test to_pascals(1.0u"MPa") ≈ 1e6
    @test to_pascals(1.0ksi) ≈ 6.894757e6 rtol=1e-5
    
    # to_newtons
    @test to_newtons(1.0u"N") ≈ 1.0
    @test to_newtons(1.0u"kN") ≈ 1000.0
    @test to_newtons(1.0kip) ≈ 4448.2216152605 rtol=1e-6
    
    # to_newton_meters
    @test to_newton_meters(1.0u"N*m") ≈ 1.0
    @test to_newton_meters(1.0u"kN*m") ≈ 1000.0
    
    # to_newtons_per_meter
    @test to_newtons_per_meter(1.0u"N/m") ≈ 1.0
    @test to_newtons_per_meter(1.0u"kN/m") ≈ 1000.0
    
    # to_kg_per_m3
    @test to_kg_per_m3(1.0u"kg/m^3") ≈ 1.0
    @test to_kg_per_m3(1.0pcf) ≈ 16.01846337 rtol=1e-5
    
    # to_m_per_s2
    @test to_m_per_s2(9.81u"m/s^2") ≈ 9.81
    @test to_m_per_s2(GRAVITY) ≈ 9.80665
end

@testset "Catalog Parsing Utilities" begin
    # asfloat
    @test asfloat(42) == 42.0
    @test asfloat(3.14) == 3.14
    @test asfloat("123.45") == 123.45
    @test_throws ArgumentError asfloat(nothing)
    
    # maybe_asfloat
    @test maybe_asfloat(42) == 42.0
    @test maybe_asfloat("123.45") == 123.45
    @test maybe_asfloat(missing) === nothing
    @test maybe_asfloat("–") === nothing  # en-dash
    @test maybe_asfloat("-") === nothing  # hyphen
    @test maybe_asfloat("—") === nothing  # em-dash
    @test maybe_asfloat("") === nothing
    @test maybe_asfloat("  ") === nothing
end

@testset "Unit Arithmetic" begin
    # Basic arithmetic should work with custom units
    @testset "Addition/Subtraction" begin
        @test 1.0kip + 1.0kip == 2.0kip
        @test 2.0ksi - 1.0ksi == 1.0ksi
        @test 100.0psf + 50.0psf == 150.0psf
    end
    
    @testset "Multiplication/Division" begin
        # Force = Pressure × Area
        @test 1.0ksi * 1.0u"inch^2" ≈ 1.0kip rtol=1e-6
        
        # Moment = Force × Length
        @test 1.0kip * 1.0u"ft" == 1.0kip*u"ft"
        
        # Linear load = Area load × Width
        @test 100.0psf * 10.0u"ft" ≈ 1000.0u"lbf/ft" rtol=1e-6
    end
    
    @testset "Unit Conversion" begin
        # Convert between compatible units
        @test uconvert(u"N", 1.0kip) ≈ 4448.2216152605u"N" rtol=1e-6
        @test uconvert(u"Pa", 1.0ksi) ≈ 6.894757e6u"Pa" rtol=1e-5
        @test uconvert(u"Pa", 1.0psf) ≈ 47.88025898u"Pa" rtol=1e-6
        @test uconvert(u"kg/m^3", 1.0pcf) ≈ 16.01846337u"kg/m^3" rtol=1e-5
    end
end
