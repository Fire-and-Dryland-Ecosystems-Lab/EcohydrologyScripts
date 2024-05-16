# global vars
site = "BigCreek"

# Ward_msr90m_basinsoilspin400y.world
# "worldfiles/Ward_msr90m_soilvegspun.world"
# reset_worlds = c("worldfiles/Ward_msr90m_world_3_thin1.world", # "worldfiles/Ward_msr90m_basinsoilspin400y_reset.world",
#                  "worldfiles/Ward_msr90m_world_3_thin2.world",
#                  "worldfiles/Ward_msr90m_world_3_thin3.world",
#                  "worldfiles/Ward_msr90m_world_3_thin4.world",
#                  "worldfiles/Ward_msr90m_world_3_thin5.world")
reset_worlds = NA
# worldfiles = c("worldfiles/Ward_msr90m_world_4_thin1.world", # Ward_msr90m_baseline_init.world", 
#                "worldfiles/Ward_msr90m_world_4_thin2.world",
#                "worldfiles/Ward_msr90m_world_4_thin3.world",
#                "worldfiles/Ward_msr90m_world_4_thin4.world",
#                "worldfiles/Ward_msr90m_world_4_thin5.world")
worldfiles = NA
# flowtables = c("flowtables/Ward_msr90m.flow",
#                "flowtables/Ward_msr90m_thin2.flow",
#                "flowtables/Ward_msr90m_thin3.flow",
#                "flowtables/Ward_msr90m_thin4.flow",
#                "flowtables/Ward_msr90m_thin5.flow")
flowtables = NA
# , "Clearcut"
# "worldfiles/Ward_msr90m_thin6_init.world"
# "flowtables/Ward_msr90m_thin6.flow"

runs = c("baseline", "thin2", "thin3", "thin4", "thin5")

runs_fullnames = c("Baseline","Prescribed Fire", "Heavy Thinning","Moderate Thinning", "Mastication")

scenario_df = data.frame(worldfiles, flowtables, runs, runs_fullnames)

# scenario_df = scenario_df[!scenario_df$runs_fullnames == "Clearcut",]

spinup_df = scenario_df
spinup_df$reset_worlds = reset_worlds
# treat_date = "1988 3 1"
# treat_date = "1996 3 1"
# dates = c("1989 9 30 24", "2020 9 30 24")
# dates = c("2005 1 1 1", "2007 9 30 24")
# dates = c("1987 9 30 24", "2004 10 1 1") #dry
# dates = c("1995 9 30 24", "2004 10 1 1") # wet
# start_output_date = "1987 10 1 1"
# start_output_date = "1995 10 1 1"

# shifting start earlier by 3 years

date_df = data.frame(climate = c("dry","wet"),
                     start = c("1984 9 30 24", "1992 9 30 24"),
                     end = c("2004 10 1 1","2012 10 1 1"),
                     start_output = c("1987 10 1 1", "1995 10 1 1"),
                     treat_date = c("1988 3 1", "1996 3 1"))

scenario_df = merge(merge(scenario_df, expand.grid(runs_fullnames = scenario_df$runs_fullnames, 
                                             climate = date_df$climate ), by = "runs_fullnames",sort = F), date_df,by = "climate",sort=F)

scenario_df = scenario_df[order(scenario_df$runs, scenario_df$climate),]

scenario_df$name = paste0("Ward_msr_", scenario_df$runs,"_", scenario_df$climate)

