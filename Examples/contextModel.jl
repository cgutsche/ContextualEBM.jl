#############################
#     Context Modeling      #
#############################

@newContext BatteryCharging
@newContext BatteryDischarging
@newContext BatteryIdle
activateContext(BatteryIdle)
BatteryChargingMode = ContextGroup(BatteryCharging, BatteryDischarging, BatteryIdle)

@newContext BatterySoCMax
@newContext BatterySoCMedium
@newContext BatterySoCMin
activateContext(BatterySoCMedium)
BatterySoC = ContextGroup(BatterySoCMax, BatterySoCMedium, BatterySoCMin)

@newContext ElectrolyzerActive
@newContext ElectrolyzerStandby
activateContext(ElectrolyzerActive)
ElectrolyzerOperation = ContextGroup(ElectrolyzerActive, ElectrolyzerStandby)

@newContext GridInputHighWithElectrolyzer
@newContext GridInputHigh
@newContext GridInputLow
@newContext GridOutput
activateContext(GridInputHigh)
GridState = ContextGroup(GridInputHighWithElectrolyzer, GridInputHigh, GridInputLow, GridOutput)

@newContext ElectrolyzerMaxPower
deactivateContext(ElectrolyzerMaxPower)
@newContext ElectrolyzerMinPower
activateContext(ElectrolyzerMinPower)

weakInclusion((!ElectrolyzerMaxPower) & BatteryCharging => BatteryIdle)
weakInclusion((!ElectrolyzerMaxPower) & (!BatterySoCMin) & (BatteryIdle) => BatteryDischarging)
weakInclusion(BatterySoCMin & BatteryDischarging => BatteryIdle)
weakInclusion(BatterySoCMax & BatteryCharging => BatteryIdle)
weakInclusion((GridOutput | GridInputLow) & ElectrolyzerActive & ElectrolyzerMinPower => ElectrolyzerStandby)
weakInclusion(GridInputHigh & ElectrolyzerStandby => ElectrolyzerActive)
weakInclusion(GridInputHighWithElectrolyzer & (ElectrolyzerMaxPower) & BatteryIdle & (!BatterySoCMax)  => BatteryCharging)

