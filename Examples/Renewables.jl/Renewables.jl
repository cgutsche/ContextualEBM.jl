using OrdinaryDiffEq, ModelingToolkit
using Symbolics
using ModelingToolkitStandardLibrary.Electrical
using ModelingToolkitStandardLibrary.Mechanical
using ModelingToolkitStandardLibrary.Mechanical.Rotational
using ModelingToolkitStandardLibrary.Blocks
using CSV, NCDatasets
using DataFrames
using Interpolations

function performanceCurve(v)
	p = 10^3 .* 2 .* [0, 0, 0, 0, 100, 200, 350, 500, 800, 1200, 1500, 1850, 2300, 2500, 2750, 2900, 3000, 3000, 3000]
	p_interp = extrapolate(interpolate(p, BSpline(Linear())), Flat())
	p_interp(max(1, v+1))
end

Symbolics.@register_symbolic performanceCurve(v)
Symbolics.@register_symbolic v_fun(t)
Symbolics.@register_symbolic i_fun(t)

@variables t
D = Differential(t)

@mtkmodel VariableResistor begin
    @extend v, i = oneport = OnePort()
    @components begin
		R = RealInput()
	end
    @equations begin
		v ~ i * R.u
    end
end

# Adapted from Modelica Standard Library
@mtkmodel Battery begin
	@extend OnePort()
	@parameters begin
		capacity = 5
		Vmax = 4.2
		Vmin = 2.5
		Idis = 0
		Ns = 1				# number of parallel connected cells
		Np = 1				# number of serial connected cells
		R_i = 0.0001
		SoC_tolerance = 10^(-9)
		SoC_0 = 0.5
	end
	@variables begin
		C(t) = SoC_0 * capacity
		SoC(t) = SoC_0, [irreducible=true]
		OCV(t) = Vmin
	end
	@components begin
		selfDischarge = Conductor(G = Np * Idis / (Ns * Vmax))
		r0 = Resistor(R = Ns * R_i / Np)
		ovc = Voltage()
		i_sensor = CurrentSensor()
		VInp = RealInput()
	end
	@equations begin
		SoC ~ C / capacity
		D(C) ~ (i_sensor.i * (SoC >= SoC_tolerance) * (i_sensor.i < 0) + i_sensor.i * (SoC <= 1-SoC_tolerance) * (i_sensor.i > 0) + 0) / Np
		OCV ~ Ns * ((Vmax - Vmin) * (0.5*(2*(SoC-0.5))^7 + 0.5) + Vmin)
		VInp.u ~ OCV
		connect(ovc.V, VInp)
		connect(p, i_sensor.p, selfDischarge.p)
		connect(ovc.p, i_sensor.n)
		connect(ovc.n, r0.p, selfDischarge.n)
		connect(n, r0.n)
	end
end

@mtkmodel BatteryWithEMS begin
	@extend OnePort()
	@parameters begin
		capacity = 5
		Vmax = 4.2
		Vmin = 2.5
		Idis = 0
		Ns = 1				# number of serial connected cells
		Np = 1				# number of parallel connected cells
		R_i = 0.0001
		SoC_tolerance = 10^(-9)
		I_charge = capacity * Ns
		SoC_max = 0.85
		SoC_min = 0.15
		SoC_0 = 0.5
		charging_0 = false
		discharging_0 = false
	end
	@variables begin
		SoC(t) = SoC_0, [irreducible=true]
		discharging(t)::Bool = discharging_0
		charging(t)::Bool = charging_0
	end
	@components begin
		battery = Battery(capacity = capacity,
		                  Vmax = Vmax,
						  Vmin = Vmin,
						  Idis = Idis,
						  Ns = Ns,
						  Np = Np,
						  R_i = R_i,
						  SoC_tolerance = SoC_tolerance)
		I_inp_grid = RealInput()
		I_inp_battery = RealInput()
		I_source_grid = Current()
		I_source_battery = Current()
		V_sensor_grid = VoltageSensor()
		V_sensor_battery = VoltageSensor()
	end
	@equations begin
		SoC ~ battery.SoC
		D(discharging) ~ 0
		D(charging) ~ 0
		I_inp_grid.u ~ battery.i_sensor.i * V_sensor_battery.v / V_sensor_grid.v
		I_inp_battery.u ~ (charging > 0.5) * I_charge - (discharging > 0.5) * I_charge
		connect(I_source_grid.I, I_inp_grid)
		connect(p, I_source_grid.n, V_sensor_grid.n)
		connect(n, I_source_grid.p, V_sensor_grid.p)
		connect(I_source_battery.I, I_inp_battery)
		connect(I_source_battery.p, battery.n, V_sensor_battery.n)
		connect(I_source_battery.n, battery.p, V_sensor_battery.p)
	end
	@continuous_events begin
		[SoC ~ SoC_max] => [charging ~ false]
		[SoC ~ SoC_min] => [discharging ~ false]
	end
	@discrete_events begin
		(SoC >= SoC_max) => [charging ~ false]
		(SoC <= SoC_min) => [discharging ~ false]
	end
end

@mtkmodel Generator begin
	@extend OnePort()
	@parameters begin
		R = 0.5 # [Ohm] armature resistance
		L = 4.5e-3 # [H] armature inductance
		k = 1 # [N.m/A] motor constant
		f = 0.01 # [N.m.s/rad] friction factor
	end
	@components begin
		R1 = Resistor(R = R)
		L1 = Inductor(L = L)
		emf = EMF(k = k)
		fixed = Fixed()
		friction = Damper(d = f)
	end
	@equations begin
		connect(emf.support, friction.flange_b, fixed.flange)
        connect(emf.flange, friction.flange_a)
		connect(p, R1.p)
		connect(R1.n, L1.p)
		connect(L1.n, emf.p)
		connect(emf.n, n)
	end
end

# Adapted from:
# Eberhart, Philip et al. (2015-09). “Open Source Library for the Simulation of Wind Power Plants”. 
# In: Proceedings of the 11th International Modelica Conference, Versailles, France,
# September 21-23, 2015. Linköping Electronic Conference Proceedings, pp. 929–936. DOI: 10.3384/ecp15118929.
@mtkmodel WindTurbine begin
    @parameters begin
        ρ=1.2
        R=90
		lambda = 7
    end
    @components begin
        tau = RealOutput()
		v_wind = RealInput()
    end
    @variables begin
        τ(t) = 0
        ϕ(t) = 0
        ω(t) = 0
        P(t) = 0
        P_Wind(t) = 0
        Vw(t) = 0
		Cp(t) = 0
    end
    @equations begin
        Vw ~ v_wind.u #v_fun(t)
		lambda * Vw ~ ω * R
		P_Wind ~ 0.5 * π * ρ * R ^ 2 * Vw ^ 3
        P ~ performanceCurve(Vw)
		Cp ~ P / P_Wind
        P ~ τ * ω
        D(ϕ) ~ ω
        tau.u ~ τ
    end
end

@mtkmodel WindPP begin
	@components begin
		windturbine = WindTurbine()
		torque = Torque()
		inertia = Inertia(J = 0.04)
		generator = Generator()
	end
	@equations begin
		connect(windturbine.tau, torque.tau)
		connect(torque.flange, inertia.flange_b)
		connect(generator.emf.flange, inertia.flange_a)
	end
end

@mtkmodel RotPowerSensor begin
    @components begin
        flange_a = Flange()
        flange_b = Flange()
        power = RealOutput()
    end
    @equations begin
        0 ~ flange_a.tau + flange_b.tau
        0 ~ flange_b.phi - flange_a.phi
        power.u ~ flange_a.tau*D(flange_a.phi)
    end
end

@mtkmodel MyInertia begin
    @parameters begin
        J, [description = "Moment of inertia"]
    end
    @components begin
        flange_a = Flange()
        flange_b = Flange()
    end
    @variables begin
        phi(t) = 0, [description = "Absolute rotation angle"]
        w(t) = 0, [description = "Absolute angular velocity"]
        a(t) = 0, [description = "Absolute angular acceleration"]
    end
    @equations begin
        phi ~ flange_a.phi
        phi ~ flange_b.phi
        D(phi) ~ w
        D(w) ~ a
        J * a ~ flange_a.tau + flange_b.tau
    end
end

@mtkmodel WindPPSimpleDC begin
	@extend OnePort()
	@structural_parameters begin
		OutputVoltage = 230
	end
	@parameters begin
		d = 100
		J = 1.7 * 10^-7
	end
	@components begin
		windpp = WindTurbine()
		v_wind = RealInput()
		P_sensor_grid = PowerSensor()
		P_turbine = RealInput()
		minusTorque = Feedback()
		P_windPP = RotPowerSensor()
		P_grid = RealInput()
		minus = Feedback()
		integ = Integrator(k = 1 / OutputVoltage / 10^-5)
		outputSource = Current()
		torque = Torque()
		fixed = Fixed()
		friction = Damper(d = d)
		inertia = MyInertia(J = J)
		gain = Gain(k = -1)
	end
	@equations begin
		connect(v_wind, windpp.v_wind)

		connect(friction.flange_b, fixed.flange)
		connect(friction.flange_a, P_windPP.flange_b)
		P_turbine.u ~ windpp.P
		connect(P_turbine, minusTorque.input1)
		connect(P_windPP.power, minusTorque.input2)
		connect(minusTorque.output, torque.tau)
		connect(torque.flange, inertia.flange_a)
		connect(P_windPP.flange_a, inertia.flange_b)

		P_grid.u ~ P_sensor_grid.power
		connect(P_windPP.power, gain.input)
		connect(gain.output, minus.input1)
		connect(P_grid, minus.input2)
		connect(minus.output, integ.input)
		connect(integ.output, outputSource.I)
		
		connect(P_sensor_grid.nv, n)
		connect(outputSource.p, n)
		connect(P_sensor_grid.pv, P_sensor_grid.pc)
		connect(P_sensor_grid.nc, outputSource.n)
		connect(P_sensor_grid.pc, p)
	end
end

# Adapted from:
# Brkic, Jovan et al. (2019-03). “Open Source PhotoVoltaics Library for Systemic Investigations”.
# In: Proceedings of the 13th International Modelica Conference, Regensburg, Germany, 
# March 4–6, 2019. Linköping Electronic Conference Proceedings. DOI: 10.3384/ecp1915741.
@mtkmodel Diode begin
	@extend v, i = oneport = OnePort()
	@parameters begin
		Bv = 18
		Ibv = 0.7
		Nbv = 0.74
		I_S = 10^(-6)
		m = 2
		VtRef = 0.02565 # k_B * T / q_electron
		R = 10^(8)
		VRef = 30.2
		IRef = 8.54
		Vt = 0.02565 # k_B * T / q_electron, here constant for T = 300 K
		IdsRef= IRef/(exp(VRef/m/VtRef))
		VNeg = m * Vt * log(Vt / VtRef)
		VBv = -m*Nbv*log(IdsRef*Nbv/Ibv)*VtRef
		VNegLin= (-VRef/m/VtRef*(Nbv*m*VtRef)) - Bv
		Ids = IRef / (exp(VRef / m / Vt) - 1)
		ns = 1
		nsModule = 1
		npModule = 1
	end
	@variables begin
		vCell(t) = 0
		vModule(t) = 0
		iModule(t) = 0
	end
	@equations begin
		vCell ~ v/ns/nsModule
		vModule ~ v/nsModule
		iModule ~ i/npModule
		i / npModule ~ (v / ns / nsModule  > VNeg ) * (Ids * (exp(v / ns / nsModule / m / Vt) - 1) + v / ns / nsModule / R) + 
		    (v / ns / nsModule  > VBv) * (v  <= VNeg) * (Ids * v / ns / nsModule / m / VtRef + v  / ns / nsModule/ R) + 
            (v / ns / nsModule  > VNegLin) * (v / ns / nsModule  <= VBv) * ((-Ibv * exp(-(v / ns / nsModule + Bv) / (Nbv * m * Vt))) + Ids * VBv / m / VtRef + v / ns / nsModule / R) +
            (v / ns / nsModule  < VNegLin) * ( Ids * v / m / Vt - Ibv * exp(VRef / m / VtRef) * (1 - (v / ns / nsModule + Bv) / (Nbv * m * Vt) - VRef / m / VtRef) + v / ns / nsModule / R)
	end
end

# Adapted from:
# Brkic, Jovan et al. (2019-03). “Open Source PhotoVoltaics Library for Systemic Investigations”.
# In: Proceedings of the 13th International Modelica Conference, Regensburg, Germany, 
# March 4–6, 2019. Linköping Electronic Conference Proceedings. DOI: 10.3384/ecp1915741.
@mtkmodel PVCell begin
	@extend OnePort()
	@parameters begin
		IscRef = 8.54
		irradianceRef = 1000
	end
	@components begin
		I_source = Current()
		I = RealInput()
		I_inp = RealInput()
		diode = Diode()
	end
	@equations begin
		I_inp.u ~ max(0, I.u / irradianceRef) * IscRef
		connect(I_source.n, p)
		connect(I_source.p, n)
		connect(I_source.I, I_inp)
		connect(I_source.p, diode.n)
		connect(I_source.n, diode.p)
	end
end

@mtkmodel PVModule begin
	@extend OnePort()
	@structural_parameters begin
		N_serial = 48
		N_parallel = 1
	end
	@components begin
		pvc_array = [PVCell() for i in 1:(N_serial*N_parallel)]
		I = RealInput()
	end
	@equations begin
		[connect(pvc_array[j].I, I) for j in 1:(N_parallel*N_serial)]...
		([[connect(pvc_array[i+(j-1)*N_serial].p, pvc_array[i+(j-1)*N_serial+1].n) for i in 1:(N_serial-1)] for j in 1:N_parallel]...)...
		[connect(pvc_array[(j-1)*N_serial+1].n, pvc_array[j*N_serial+1].n) for j in 1:(N_parallel-1)]...
		[connect(pvc_array[j*N_serial].p, pvc_array[(j+1)*N_serial].p) for j in 1:(N_parallel-1)]...
		connect(pvc_array[1].n, n)
		connect(pvc_array[N_serial].p, p)
	end
end

@mtkmodel PVModuleSimple begin
	@extend OnePort()
	@structural_parameters begin
		nsModule = 1
		npModule = 1
	end
	@parameters begin
		IscRef = 8.54
		irradianceRef = 1000
		ns = 1
	end
	@components begin
		I_source = Current()
		I = RealInput()
		I_inp = RealInput()
		diode = Diode(ns = ns, nsModule = nsModule, npModule = npModule)
	end
	@equations begin
		I_inp.u ~ max(0, I.u / irradianceRef) * IscRef
		connect(I_source.n, p)
		connect(I_source.p, n)
		connect(I_source.I, I_inp)
		connect(I_source.p, diode.n)
		connect(I_source.n, diode.p)
	end
end

@mtkmodel PVPark_simplified begin
	@extend OnePort()
	@structural_parameters begin
		P_peak_installed = 400
		P_peak_ref = 1000
		OutputVoltage = 230
	end
	@components begin
		P_sensor_grid = PowerSensor()
		P_PV =  RealInput()
		P_grid = RealInput()
		minus = Feedback()
		integ = Integrator(k = 1 / OutputVoltage / 10^-5)
		outputSource = Current()
		gain = Gain(k = -1)
		iridiation_input = RealInput()
	end
	@variables begin
		iridiation(t) = 0
		P(t) = 0
	end
	@equations begin
		iridiation ~ iridiation_input.u #i_fun(t) # D(iridiation) ~ (t > 25) * 20 - (t > 75) * 20 # iridiation ~ i_fun(t)
		P ~ max(0, min(1, iridiation/P_peak_ref))*P_peak_installed
		
		P_PV.u ~ P
		P_grid.u ~ P_sensor_grid.power
		connect(P_PV, gain.input)
		connect(gain.output, minus.input1)
		connect(P_grid, minus.input2)
		connect(minus.output, integ.input)
		connect(integ.output, outputSource.I)
		
		connect(P_sensor_grid.nv, n)
		connect(outputSource.p, n)
		connect(P_sensor_grid.pv, P_sensor_grid.pc)
		connect(P_sensor_grid.nc, outputSource.n)
		connect(P_sensor_grid.pc, p)
	end
end

# Adapted from:
# Migoni G, et al., "Efficient simulation of Hybrid Renewable Energy Systems", 
# International Journal of Hydrogen Energy (2016), http://dx.doi.org/10.1016/j.ijhydene.2016.06.019
@mtkmodel PEMElectrolyzerStack begin
	@extend v, i = oneport = OnePort()
	@parameters begin
		F = 96485 # [C mol^-1]
		v_rev = 1.229 # [V]
		r1 = 7.331e-5 # [ohm m^2]
		r2 = -1.107e-7 # [ohm m^2 C-1]
		r3 = 0
		s1 = 1.586e-1 # [C]
		s2 = 1.378e-3 # [V C-1]
		s3 = -1.606e-5 # [V C-2]
		t1 = 1.599e-2 # [m^2 A^-1]
		t2 = -1.302 # [m^2 A-1 C-1]
		t3 = 4.213e2 # [m^2 A-1 C-2]
		A = 0.25 # [m2]
		nc = 1 #number of serial cells
		T = 40 # 40
	end
	@components begin
		y = RealOutput()
	end
	@variables begin
		η(t) = 0	# Faraday efficiency
		v_cell(t) = 0
		v_res(t) = 0
		v_act(t) = 0
		P(t) = 0
		fH2(t) = 0
	end
	@equations begin
		η ~ (-0.00075 * T + 1) * ((i / A)^2)/( (2.5 * T + 50) + (i / A)^2)
		P ~ i * v
		v ~ v_cell * nc
		v_cell ~ v_rev + v_res + v_act
		v_res ~ ((r1 + r2 * T) * i / A)
		v_act ~ (s1 + s2 * T + s3 * T^2) * log10(((t1 + t2 / T + t3 / T^2) * i / A) + 1)
		fH2 ~ η * P / (v * 2 * F) * 3600
		y.u ~ fH2
	end
end

@mtkmodel PEMElectrolyzer begin
	@extend OnePort()
	@structural_parameters begin
		n_stacks = 1
		n_cells_per_stack = 1 #number of serial cells
		T = 80 # 40
	end
	@components begin
		stacks = [PEMElectrolyzerStack(nc = n_cells_per_stack) for i in 1:n_stacks]
		H2output = RealInput()
	end
	@variables begin
		η(t) = 0	# Faraday efficiency
		v_cell(t) = 0
		v_res(t) = 0
		v_act(t) = 0
		P(t) = 0
		fH2(t) = 0
	end
	@equations begin
		[connect(stacks[i].n, stacks[i+1].n) for i in 1:(n_stacks-1)]...
		[connect(stacks[i].p, stacks[i+1].p) for i in 1:(n_stacks-1)]...
		connect(n, stacks[1].n)
		connect(p, stacks[1].p)
		H2output.u ~ sum([stacks[i].y.u for i in 1:n_stacks])
	end
end

# Adapted from:
# Migoni G, et al., "Efficient simulation of Hybrid Renewable Energy Systems", 
# International Journal of Hydrogen Energy (2016), http://dx.doi.org/10.1016/j.ijhydene.2016.06.019
@mtkmodel FuelCellStack begin
	@extend v, i = oneport = OnePort()
    @parameters begin
        p_H2 = 0.3           # Pressure of hydrogen [atm]
        p_O2 = 1             # Pressure of oxygen [atm]
        Nc = 47              # Number of series cells
        R = 8.3143           # Gas constant [J/(mol*K)]
        F = 96487            # Faraday constant [C/mol]
        ke = 8.5e-4          # Constant [V/K]
        deltaG = -237153.66  # Gibbs free energy per mole of reaction [J/mol]
        ne = 2               # Number of electrons in reaction
        nu_0 = 26.5230       # Voltage constant [V]
        a = 8.9224e-2        # Activation energy constant
        I_limit = 75         # Fuel cell current limit [A]
        T = 298              # Temperature [K]
        MH2 = 2.016e-3       # Molar mass of hydrogen [kg/mol]
    end
    @variables begin
        HF2_FC(t) = 0         # Hydrogen flow [kg/s]
        EO_cell_std(t) = 0    # Standard reference potential [V]
        EO_cell(t) = 0        # Reference potential [V]
        E_cell(t) = 0         # Nernst Equation (ideal voltage) [V]
        EStack(t) = 0         # Total stack voltage [V]
        Vohm(t) = 0           # Ohmic voltage drop [V]
        Vconc(t) = 0          # Concentration voltage drop [V]
        Vact(t) = 0           # Activation voltage drop [V]
        React0(t) = 0         # Reaction constant term 0
        React1(t) = 0         # Reaction constant term 1
        React2(t) = 0         # Reaction constant term 2
    end
    @equations begin
		EO_cell_std ~ -deltaG / (ne * F)
        EO_cell ~ EO_cell_std + ke * (T - 298)
        E_cell ~ EO_cell + (R * T / (2 * F)) * log(p_H2 * sqrt(p_O2))
        EStack ~ E_cell * Nc
        React0 ~ -1.0526
        React1 ~ (6.945e-11) * i^6 - (1.7272e-8) * i^5 + (1.7772e-6) * i^4 -
                 (9.8133e-5) * i^3 + (3.1430e-3) * i^2 - (3.5320e-2) * i
        React2 ~ (1.389e-3) * (T - 298)
        Vact ~ -nu_0 - (T - 298) * (a * (React0 + React1 + React2) + i)
        Vohm ~ i * (1.7941 - (2.308e-2) * i - (2.006e-3) * (T - 298))
        Vconc ~ -(R * T / (ne * F)) * log(1 - (i / I_limit))
        v ~ EStack - Vact - Vohm - Vconc
        HF2_FC ~ MH2 * Nc * i / (2 * F) * 3600
	end
end

# Adapted from:
# Migoni G, et al., "Efficient simulation of Hybrid Renewable Energy Systems", 
# International Journal of Hydrogen Energy (2016), http://dx.doi.org/10.1016/j.ijhydene.2016.06.019
@mtkmodel FuelCell begin
	@extend v, i = oneport = OnePort()
    @parameters begin
        p_H2 = 0.3           # Pressure of hydrogen [atm]
        p_O2 = 1             # Pressure of oxygen [atm]
        Nc = 47              # Number of series cells
		Nstacks = 10         # Number of series cells
        R = 8.3143           # Gas constant [J/(mol*K)]
        F = 96487            # Faraday constant [C/mol]
        ke = 8.5e-4          # Constant [V/K]
        deltaG = -237153.66  # Gibbs free energy per mole of reaction [J/mol]
        ne = 2               # Number of electrons in reaction
        nu_0 = 26.5230       # Voltage constant [V]
        a = 8.9224e-2        # Activation energy constant
        I_limit = 75         # Fuel cell current limit [A]
        T = 298              # Temperature [K]
        MH2 = 2.016e-3       # Molar mass of hydrogen [kg/mol]
    end
    @variables begin
        HF2_FC(t) = 0         # Hydrogen flow [kg/s]
        EO_cell_std(t) = 0    # Standard reference potential [V]
        EO_cell(t) = 0        # Reference potential [V]
        E_cell(t) = 0         # Nernst Equation (ideal voltage) [V]
        EStack(t) = 0         # Total stack voltage [V]
        Vohm(t) = 0           # Ohmic voltage drop [V]
        Vconc(t) = 0          # Concentration voltage drop [V]
        Vact(t) = 0           # Activation voltage drop [V]
        React0(t) = 0         # Reaction constant term 0
        React1(t) = 0         # Reaction constant term 1
        React2(t) = 0         # Reaction constant term 2
		IStack(t) = 0         # Reaction constant term 2
	end
	@components begin
		H2input = RealInput()
	end
    @equations begin
		EO_cell_std ~ -deltaG / (ne * F)
        EO_cell ~ EO_cell_std + ke * (T - 298)
        E_cell ~ EO_cell + (R * T / (2 * F)) * log(p_H2 * sqrt(p_O2))
        EStack ~ E_cell * Nc
        React0 ~ -1.0526
        React1 ~ (6.945e-11) * IStack^6 - (1.7272e-8) * IStack^5 + (1.7772e-6) * IStack^4 -
                 (9.8133e-5) * IStack^3 + (3.1430e-3) * IStack^2 - (3.5320e-2) * IStack
        React2 ~ (1.389e-3) * (T - 298)
        Vact ~ nu_0 / 47 * Nc - (T - 298) * (a * (React0 + React1 + React2) + IStack)
        Vohm ~ IStack * (1.7941 - (2.308e-2) * IStack - (2.006e-3) * (T - 298)) / 47 * Nc 
        Vconc ~ -(R * T / (ne * F)) * log(1 - (IStack / I_limit))
        v ~ EStack - Vact - Vohm - Vconc
		IStack ~ i / Nstacks
        HF2_FC ~ MH2 * Nc * i / (2 * F) * 3600
		HF2_FC ~ H2input.u
	end
end

@mtkmodel FuelCellSimple begin
	@extend v, i = oneport = OnePort()
    @parameters begin
        Nc = 150               # Number of series cells
		Nstacks = 10
		V_cell = 0.7
		MH2 = 2.016e-3       # Molar mass of hydrogen [kg/mol]
		F = 96487            # Faraday constant [C/mol]
    end
    @variables begin
		I_cell(t) = 0
		F_H2_FC(t) = 0
    end
	@components begin
		H2input = RealInput()
	end
    @equations begin
        v ~ V_cell * Nc
		I_cell ~ i / Nstacks 
        F_H2_FC ~ MH2 * Nc * i / (2 * F) * 3600
		H2input.u ~ F_H2_FC
	end
end

@mtkmodel FuelCellWithDCDC begin
	@extend OnePort()
    @parameters begin
        Nc = 150               # Number of series cells
		Nstacks = 10
		V_cell = 0.7
		MH2 = 2.016e-3       # Molar mass of hydrogen [kg/mol]
		F = 96487            # Faraday constant [C/mol]
		P = 1000000
    end
	@components begin
		fc = FuelCellSimple()
		H2input = RealInput()
		outputCurrent = Current()
		VSensor = VoltageSensor()
		R = Resistor(R = (V_cell*Nc)^2/P)
		IInp = RealInput()
	end
    @equations begin
        connect(fc.p, R.n)
	    connect(fc.n, R.p)

		connect(outputCurrent.p, VSensor.p, n)
		connect(outputCurrent.n, VSensor.n, p)

		IInp.u ~ fc.i * fc.v / VSensor.v
		connect(outputCurrent.I, IInp)
	end
end

@mtkmodel DCDC_Converter begin
	@components begin
		in = OnePort()
		out = OnePort()
		outputVoltage = RealInput()
	end
	@parameters begin
		η(t) = 1.0
	end
	@equations begin
		out.v ~ outputVoltage.u
		out.i ~ - η * in.v * in.i / out.v
	end
end

# Adapted from:
# Migoni G, et al., "Efficient simulation of Hybrid Renewable Energy Systems", 
# International Journal of Hydrogen Energy (2016), http://dx.doi.org/10.1016/j.ijhydene.2016.06.019
@mtkmodel H2Tank begin
    @parameters begin
        R = 8.314  # Universal gas constant [J/(mol*K)]
        a = 0.0247 # van der Waals coefficient for H2 [Pa·m6/mol2]
        b = 0.0266 # van der Waals coefficient for H2 [m3/mol]
        T = 300    # Constant temperature assumption [K]
        M_H2 = 0.002  # Molar mass of H2 [kg/mol]
        V_tank = 0.1 # Volume of the single storage tank [m3]
    end

    @variables begin
        H2Mass(t)      # Hydrogen mass in the tank [kg]
        H2Moles(t)   # Number of hydrogen moles in the tank
        Pressure(t)  # Pressure inside the tank [Pa]
        SoC(t)  # State of charge (SoC) of the tank
    end

    @components begin
        inflow = RealInput()   # Hydrogen input to the tank [kg/s]
        outflow = RealInput()  # Hydrogen output from the tank [kg/s]
        integrator = Integrator()  # Integrator for mass balance
    end

    @equations begin
        H2Moles ~ H2Mass / M_H2

        # van der Waals Equation
        Pressure ~ (R * T * H2Moles) / (V_tank - H2Moles * b) - (a * H2Moles^2) / V_tank^2

        # State of Charge (SoC)
        SoC ~ Pressure / 8e6  # Normalize with max pressure of 80 bar (8e6 Pa)

        # Mass balance: inflow - outflow = change in hydrogen mass
        connect(inflow, integrator.input)
        connect(outflow, integrator.output)
        H2Mass ~ integrator.y
    end
end