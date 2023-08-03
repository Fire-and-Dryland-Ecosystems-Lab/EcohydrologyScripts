# output aliase
#soils spinup
output_vars_soils = c("patch.soil_cs.totalc", "patch.soil_cs.DOC", "patch.soil_ns.totaln", "patch.soil_ns.DON", "stratum.cs.net_psn", 
                      "patch.sat_deficit", "patch.lai", "patch.totalc")

output_vars_soils_debug = c("patch.soil_cs.totalc", "patch.soil_cs.DOC", "patch.soil_ns.totaln", "patch.soil_ns.DON", "stratum.cs.net_psn", "patch.sat_deficit", "patch.lai", 
                            "patch.totalc", "patch.evaporation","patch.evaporation_surf", "patch.transpiration_unsat_zone","patch.transpiration_sat_zone",  "stratum.cs.net_psn",
                            "patch.snowpack.water_equivalent_depth", "patch.rootzone.depth", "patch.gw_drainage", "patch.rz_storage", "patch.unsat_storage")

# streamflow calibration
output_vars_streamflowcal = c("patch.evaporation", "patch.transpiration_unsat_zone", "patch.transpiration_sat_zone", "patch.streamflow",
                              "patch.snowpack.water_equivalent_depth", "patch.sat_deficit", "patch.lai", "stratum.epv.height", "patch.rootzone.depth" )

# water balance +
output_vars_waterbal = c("patch.evaporation","patch.evaporation_surf", "patch.exfiltration_sat_zone", "patch.exfiltration_unsat_zone",
                         "patch.transpiration_unsat_zone","patch.transpiration_sat_zone", "patch.streamflow",
                         "patch.totalc", "stratum.cs.net_psn", "stratum.cdf.psn_to_cpool", "patch.snowpack.water_equivalent_depth", "hill.gw.Qout",
                         "patch.rootzone.depth", "patch.sat_deficit", "patch.gw_drainage", "patch.lai","patch.rz_storage",
                         "patch.rz_transfer", "patch.unsat_storage", "patch.unsat_transfer", "patch.detention_store",
                         "stratum.Kstar_potential_both", "patch.total_water_in", "patch.litter.rain_stored", "hill.gw.storage","patch.canopy_rain_stored",
                         "patch.canopy_snow_stored")

output_vars_minimal = c("patch.soil_cs.totalc", "patch.soil_ns.totaln", "stratum.cs.net_psn", 
                        "patch.sat_deficit", "patch.lai", "patch.totalc", "patch.evaporation", "patch.streamflow", "patch.total_water_in")

output_vars_interest = c("stratum.cs.net_psn", "patch.sat_deficit", "patch.lai", "patch.totalc", "patch.evaporation", "patch.streamflow")

output_cpools = c("stratum.cs.totalc", "stratum.cs.leafc", "stratum.cs.cpool", "stratum.cs.live_stemc", "stratum.cs.dead_stemc", "stratum.cs.live_crootc", "stratum.cs.dead_crootc", "stratum.cs.frootc")
