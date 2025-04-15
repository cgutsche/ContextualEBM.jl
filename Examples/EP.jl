using Pkg
Pkg.activate("../../../ContextualJulia/Test6")
using Plots, Parameters, Symbolics
pyplot()
using BenchmarkTools
using OrdinaryDiffEq, ModelingToolkit
using CSV, NCDatasets
using Interpolations

using Contexts
include("../CVSS.jl")
include("Renewables.jl/Renewables.jl")


#############################
#     Context Modeling      #
#############################

include("contextModel.jl")

#############################
#     Data preparation      #
#############################

file = open("solarData.txt")
data = readlines(file)
radData = parse.(Float64, data)

file = open("windData.txt")
data = readlines(file)
windData = parse.(Float64, data)

v_interp = interpolate(windData, BSpline(Quadratic(Reflect(OnCell()))))
i_interp = interpolate(radData, BSpline(Linear()))

v_fun(t) = v_interp(Float64(t)+1)
i_fun(t) = i_interp(Float64(t)+1)

Symbolics.@register_symbolic v_fun(t)
Symbolics.@register_symbolic i_fun(t)

#############################
#            COP            #
#############################

@context BatteryCharging function batteryMode(integ, u)
    integ.u[u.charging] = true
    integ.u[u.discharging] = false
end

@context BatteryIdle function batteryMode(integ, u)
    integ.u[u.charging] = false
    integ.u[u.discharging] = false
end

@context BatteryDischarging function batteryMode(integ, u)
    integ.u[u.charging] = false
    integ.u[u.discharging] = true
end

#############################
#         Callbacks         #
#############################

function changeGrid!(integ, u, p, ctx)
    oldEloContext = ElectrolyzerOperation()
    oldBatteryContext = BatteryChargingMode()
    if isActive(ElectrolyzerActive)
        if (integ.u[u.P_grid] <= 0)
            activateContext(GridOutput) 
        elseif (integ.u[u.P_grid] <= 100)
            activateContext(GridInputLow) 
        elseif (integ.u[u.P_grid] <= 3000000) & isActive(GridInputLow | GridOutput)
            activateContext(GridInputHigh) 
        elseif (integ.u[u.P_grid] <= 10000) & isActive(GridInputHighWithElectrolyzer)
            activateContext(GridInputHigh) 
        elseif (integ.u[u.P_grid] >= 3000000)
            if isActive(GridInputHigh & (ElectrolyzerMaxPower) & (BatteryIdle)) & (integ.t > 0.01)
                activateContext(GridInputHighWithElectrolyzer)
            end
            if isActive(GridInputHigh & (ElectrolyzerMaxPower) & (BatteryDischarging)) & (integ.t > 0.01)
                activateContext(BatteryIdle)
            end
        end
    else
        activateContext(GridInputHigh)
    end
    if oldBatteryContext != BatteryChargingMode()
        @context BatteryChargingMode() batteryMode(integ, u)
    end
    if oldEloContext != ElectrolyzerOperation()
        restart!(integ)
    end
end

function changeBatterySoCState!(integ, u, p, ctx)
    oldEloContext = ElectrolyzerOperation()
    oldContext = BatteryChargingMode()
    if (integ.u[u.SoC] < 0.16) & isActive(BatterySoCMedium)
        activateContext(BatterySoCMin)
    elseif (integ.u[u.SoC] > 0.84) & isActive(BatterySoCMedium)
        activateContext(BatterySoCMax)
    elseif (integ.u[u.SoC] > 0.16) | (integ.u[u.SoC] < 0.84)
        activateContext(BatterySoCMedium)
    end
    if oldContext != BatteryChargingMode()
        @context BatteryChargingMode() batteryMode(integ, u)
    end
    if oldEloContext != ElectrolyzerOperation()
        restart!(integ)
    end
end

function electrolyzerMax!(integ, u, p, ctx)
    oldContext = ElectrolyzerOperation()
    if !(isActive(ElectrolyzerMaxPower)) & (integ.t > 0.01)
        activateContext(ElectrolyzerMaxPower)
        @context BatteryChargingMode() batteryMode(integ, u)
        if oldContext != ElectrolyzerOperation()
            restart!(integ)
        end
    end
end

function electrolyzerNotMax!(integ, u, p, ctx)
    oldContext = ElectrolyzerOperation()
    if isActive(ElectrolyzerMaxPower) & (integ.t > 0.01)
        deactivateContext(ElectrolyzerMaxPower)
        @context BatteryChargingMode() batteryMode(integ, u)
        if oldContext != ElectrolyzerOperation()
            restart!(integ)
        end
    end
end

function electrolyzerMin!(integ, u, p, ctx)
    oldContext = ElectrolyzerOperation()
    if !(isActive(ElectrolyzerMinPower))# & (integ.t > 0.01)
        activateContext(ElectrolyzerMinPower)
        @context BatteryChargingMode() batteryMode(integ, u)
        if oldContext != ElectrolyzerOperation()
            restart!(integ)
        end
    end
end

function electrolyzerNotMin!(integ, u, p, ctx)
    oldContext = ElectrolyzerOperation()
    if isActive(ElectrolyzerMinPower) & (integ.t > 0.01)
        deactivateContext(ElectrolyzerMinPower)
        @context BatteryChargingMode() batteryMode(integ, u)
        if oldContext != ElectrolyzerOperation()
            restart!(integ)
        end
    end
end

#############################
#         EBM Models        #
#############################

@mtkmodel EP_ElectrolyzerActive begin
    @structural_parameters begin
            N_W = 6
            P_inst = 10*10^6
        end
    @components begin
            elo_v_inp = RealInput()
            dcdc_electrolyzer = DCDC_Converter()
            i_inp = RealInput()
            pv1 =  PVPark_simplified(P_peak_installed = P_inst)
            v_wind = RealInput()
            H2 = Integrator()
            wp = [WindPPSimpleDC() for i in 1:N_W]
            battery = BatteryWithEMS(Ns = 125, Np = 500, capacity = 250, Vmax = 3.7, Vmin = 2.9, I_charge=6250)
            ground = Ground()
            V1 = RealInput()
            source = Voltage()
            bus = Pin()
            i_sensor = CurrentSensor()
            electrolyzer = PEMElectrolyzer(n_stacks = 5, n_cells_per_stack = 300)
        end
    @variables begin
        P_grid(t) = 0, [irreducible=true]
        P_wp(t) = 0, [irreducible=true]
        P_pv(t) = 0, [irreducible=true]
        P_elo(t) = 0, [irreducible=true]
        P_bat(t) = 0, [irreducible=true]
        windSpeed(t) = 0, [irreducible=true]
        V_elo(t) = 300*1.8
    end
    @equations begin
        windSpeed ~ v_fun(t)
        i_inp.u ~ i_fun(t)
        v_wind.u ~ windSpeed

        V1.u ~ 230
        connect(source.V, V1)
        connect(i_sensor.p, source.n)
        connect(i_sensor.n, bus)
        connect(source.p, ground.g)
        P_grid ~ 230 * i_sensor.i
        connect(dcdc_electrolyzer.outputVoltage, elo_v_inp)
        connect(dcdc_electrolyzer.in.n, bus)
        connect(dcdc_electrolyzer.in.p, ground.g)
        connect(electrolyzer.n, dcdc_electrolyzer.out.n)
        connect(dcdc_electrolyzer.out.p, ground.g)
        connect(electrolyzer.p, ground.g)
        P_elo ~ electrolyzer.v * electrolyzer.i
        connect(pv1.iridiation_input, i_inp)
        connect(pv1.p, bus)
        connect(pv1.n, ground.g)
        P_pv ~ pv1.v * pv1.i
        [connect(wp[i].v_wind, v_wind) for i in 1:N_W]...
        [connect(wp[i].p, bus) for i in 1:N_W]...
        [connect(wp[i].n, ground.g) for i in 1:N_W]...
        P_wp ~ sum([wp[i].v * wp[i].i  for i in 1:N_W])
        connect(battery.n, bus)
        connect(battery.p, ground.g)
        P_bat ~ battery.i * battery.v

        D(V_elo) ~ 200*(sign(P_grid - 100000)) * ((V_elo < 300*2.6) * (sign(P_grid - 200000) == 1) + (V_elo > 300*1.8) * (sign(P_grid - 50000) == -1))
        elo_v_inp.u ~ V_elo 
        connect(H2.input, electrolyzer.H2output)
    end
    @continuous_events begin
        ([battery.SoC ~ 0.155] => (changeBatterySoCState!, [battery.SoC => :SoC, battery.charging => :charging, battery.discharging => :discharging], [], []))
        ([battery.SoC ~ 0.845] => (changeBatterySoCState!, [battery.SoC => :SoC, battery.charging => :charging, battery.discharging => :discharging], [], []))
        ([battery.SoC ~ 0.165] => (changeBatterySoCState!, [battery.SoC => :SoC, battery.charging => :charging, battery.discharging => :discharging], [], []))
        ([battery.SoC ~ 0.835] => (changeBatterySoCState!, [battery.SoC => :SoC, battery.charging => :charging, battery.discharging => :discharging], [], []))
        ([P_grid  ~ 100] => (changeGrid!, [P_grid, battery.charging => :charging, battery.discharging => :discharging], [], []))
        ([P_grid  ~ 10000] => (changeGrid!, [P_grid, battery.charging => :charging, battery.discharging => :discharging], [], []))
        ([P_grid  ~ 3000000] => (changeGrid!, [P_grid, battery.charging => :charging, battery.discharging => :discharging], [], []))
        ([V_elo ~ 300*2.6-5] => (electrolyzerMax!, [battery.charging => :charging, battery.discharging => :discharging], [], []))
        ([V_elo ~ 300*2.6-10] => (electrolyzerNotMax!, [battery.charging => :charging, battery.discharging => :discharging], [], []))
        ([V_elo ~ 300*1.8+5] => (electrolyzerMin!, [battery.charging => :charging, battery.discharging => :discharging], [], []))
        ([V_elo ~ 300*1.8+10] => (electrolyzerNotMin!, [battery.charging => :charging, battery.discharging => :discharging], [], []))
    end
    @discrete_events begin
        ((P_grid < 0) => (changeGrid!, [P_grid, battery.charging => :charging, battery.discharging => :discharging], [], []))
        ((P_grid > 3000000) => (changeGrid!, [P_grid, battery.charging => :charging, battery.discharging => :discharging], [], []))
    end
end

@mtkmodel EP_ElectrolyzerStandby begin
    @structural_parameters begin
            N_W = 6
            P_inst = 10*10^6
        end
    @components begin
            elo_v_inp = RealInput()
            dcdc_electrolyzer = DCDC_Converter()
            i_inp = RealInput()
            pv1 =  PVPark_simplified(P_peak_installed = P_inst)
            v_wind = RealInput()
            H2 = Integrator()
            wp = [WindPPSimpleDC() for i in 1:N_W]
            battery = BatteryWithEMS(Ns = 125, Np = 500, capacity = 250, Vmax = 3.7, Vmin = 2.9, I_charge=6250)
            ground = Ground()
            V1 = RealInput()
            source = Voltage()
            bus = Pin()
            i_sensor = CurrentSensor()
            electrolyzer = Resistor(R = 0.66125)
            H2_inp = RealInput()
        end
    @variables begin
        P_grid(t) = 0, [irreducible=true]
        P_wp(t) = 0, [irreducible=true]
        P_pv(t) = 0, [irreducible=true]
        P_elo(t) = 0, [irreducible=true]
        P_bat(t) = 0, [irreducible=true]
        windSpeed(t) = 0, [irreducible=true]
        V_elo(t) = 300*1.8
    end
    @equations begin
        windSpeed ~ v_fun(t)
        i_inp.u ~ i_fun(t)
        v_wind.u ~ windSpeed

        V1.u ~ 230
        connect(source.V, V1)
        connect(i_sensor.p, source.n)
        connect(i_sensor.n, bus)
        connect(source.p, ground.g)
        P_grid ~ 230 * i_sensor.i
        connect(dcdc_electrolyzer.outputVoltage, elo_v_inp)
        connect(dcdc_electrolyzer.in.n, bus)
        connect(dcdc_electrolyzer.in.p, ground.g)
        connect(electrolyzer.n, dcdc_electrolyzer.out.n)
        connect(dcdc_electrolyzer.out.p, ground.g)
        connect(electrolyzer.p, ground.g)
        P_elo ~ electrolyzer.v * electrolyzer.i
        connect(pv1.iridiation_input, i_inp)
        connect(pv1.p, bus)
        connect(pv1.n, ground.g)
        P_pv ~ pv1.v * pv1.i
        [connect(wp[i].v_wind, v_wind) for i in 1:N_W]...
        [connect(wp[i].p, bus) for i in 1:N_W]...
        [connect(wp[i].n, ground.g) for i in 1:N_W]...
        P_wp ~ sum([wp[i].v * wp[i].i  for i in 1:N_W])
        connect(battery.n, bus)
        connect(battery.p, ground.g)
        P_bat ~ battery.i * battery.v 
        V_elo ~ 300*1.8
        elo_v_inp.u ~ 230
        H2_inp.u ~ 0
        connect(H2.input, H2_inp)
    end
    @continuous_events begin
        ([battery.SoC ~ 0.155] => (changeBatterySoCState!, [battery.SoC => :SoC, battery.charging => :charging, battery.discharging => :discharging], [], []))
        ([battery.SoC ~ 0.845] => (changeBatterySoCState!, [battery.SoC => :SoC, battery.charging => :charging, battery.discharging => :discharging], [], []))
        ([battery.SoC ~ 0.165] => (changeBatterySoCState!, [battery.SoC => :SoC, battery.charging => :charging, battery.discharging => :discharging], [], []))
        ([battery.SoC ~ 0.835] => (changeBatterySoCState!, [battery.SoC => :SoC, battery.charging => :charging, battery.discharging => :discharging], [], []))
        ([P_grid ~ 1000000] => (changeGrid!, [P_grid, battery.charging => :charging, battery.discharging => :discharging], [], []))
    end
    @discrete_events begin
        ((P_grid > 1000000) => (changeGrid!, [P_grid, battery.charging => :charging, battery.discharging => :discharging], [], []))
    end
end

#############################
#         ODE Systems       #
#############################

@mtkbuild sys_EP_ElectrolyzerActive = EP_ElectrolyzerActive()
@mtkbuild sys_EP_ElectrolyzerStandby = EP_ElectrolyzerStandby()

#############################
#        CVSS Problem       #
#############################

prob = CVSSProblem(Dict(ElectrolyzerActive => sys_EP_ElectrolyzerActive,
                        ElectrolyzerStandby=> sys_EP_ElectrolyzerStandby),
                    [],
                    (0, 1000))

#############################
#          Solving          #
#############################

println("Simulation solving")

bt = @elapsed begin
        global sol = solve(prob, Rodas5(autodiff=false); dtmax=0.005)
    end
println(bt)

#############################
#         Plotting          #
#############################

p_1 = plot(sol[t], sol[sys_EP_ElectrolyzerActive.P_wp]./1000000, label="Wind Turbines", ylabel = "Power [MW]", legend = :topright, xlabelfontsize=16, ylabelfontsize=16, xtickfontsize=16, ytickfontsize=16,legendfontsize=16)
p_1 = xlabel!("time [h]")
p_1 = plot!(sol[t], sol[sys_EP_ElectrolyzerActive.P_pv]./1000000, label = "PV", color=:red, legend = :topright,legendfontsize=16)
p_1 = plot!(sol[t], sol[sys_EP_ElectrolyzerActive.P_elo]./1000000, label = "Electrolzyer", color=:green, legend = :topright,legendfontsize=16)
p_1 = plot!(sol[t], sol[sys_EP_ElectrolyzerActive.i_sensor.i] .* 230 ./1000000, label="Grid", color=:purple, legend = :topright,legendfontsize=16)
p_1 = plot!(sol[t], sol[sys_EP_ElectrolyzerActive.P_bat]./1000000, label = "Battery", color=:orange, legend = :topright,legendfontsize=16)

display(p_1)
readline()

p_2 = plot(sol[t], sol[sys_EP_ElectrolyzerActive.H2.y]./1000, ylabel = "Hydrogen mass [t]", legend = :none, xlabelfontsize=16, ylabelfontsize=16)
p_2 = xlabel!("time [h]")

display(p_2)
readline()

p_3 = plot(sol[t], sol[sys_EP_ElectrolyzerActive.battery.SoC], ylabel = "battery SoC", legend = :none, xlabelfontsize=16, ylabelfontsize=16)
p_3 = xlabel!("time [h]")

display(p_3)
readline()
