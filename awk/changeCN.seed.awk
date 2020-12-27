{a = 0;}
($2 == "gw.storage") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "gw.NO3") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "lna") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "Ksat_vertical") {printf("%f  %s\n" ,1.0,$2); a=1;}
($2 == "mpar") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "rz_storage") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "unsat_storage") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "sat_deficit") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "snowpack.water_equivalent_depth") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "snowpack.water_depth") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "snowpack.T") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "snowpack.surface_age") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "snowpack.energy_deficit") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "litter.cover_fraction") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "litter.rain_stored") {printf("%f  %s\n" ,0.0000,$2); a=1;}
($2 == "litter_cs.litr1c") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "litter_ns.litr1n") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "litter_cs.litr2c") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "litter_cs.litr3c") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "litter_cs.litr4c") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "soil_cs.soil1c") {printf("%f	 %s\n",-9999,$2); a=1;}
($2 == "soil_ns.sminn") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "soil_ns.nitrate") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "soil_cs.soil2c") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "soil_cs.soil3c") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "soil_cs.soil4c") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "cover_fraction") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "gap_fraction") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "rootxone.depth") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "snow_stored") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "rain_stored") {printf("%f  %s\n" ,-9999,$2); a=1;}

($2 == "cs.cpool") {printf("%f	%s\n",0.004,$2); a=1;}
($2 == "cs.leafc") {printf("%f	%s\n",$1*0.02,$2); a=1;}
($2 == "cs.dead_leafc") {printf("%f	%s\n",$1*0.02,$2); a=1;}
($2 == "cs.leafc_store") {printf("%f	%s\n",$1*0.02,$2); a=1;}
($2 == "cs.leafc_transfer") {printf("%f	%s\n",$1*0.02,$2); a=1;}
($2 == "cs.live_stemc") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "cs.live_stemc_store") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "cs.livestemc_store") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "cs.live_stemc_transfer") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "cs.dead_stemc") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "cs.deadstemc_store") {printf("%f     %s\n",$1*0.00,$2); a=1;}
($2 == "cs.deadstemc_transfer") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "cs.live_crootc") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "cs.livecrootc_store") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "cs.livecrootc_transfer") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "cs.dead_crootc") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "cs.deadcrootc_store") {printf("%f  %s\n" ,-9999,$2); a=1;} 
($2 == "cs.deadcrootc_transfer") {printf("%f  %s\n" ,-9999,$2); a=1;} 
($2 == "cs.frootc") {printf("%f	%s\n",$1*0.02,$2); a=1;}
($2 == "cs.frootc_store") {printf("%f	%s\n",$1*0.02,$2); a=1;}
($2 == "cs.frootc_transfer") {printf("%f	%s\n",$1*0.02,$2); a=1;}
($2 == "cs.cwdc") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "epv.prev_leafcalloc") {printf("%f  %s\n" ,-9999,$2); a=1;}

($2 == "ns.leafn") {printf("%f	%s\n",$1*0.015,$2); a=1;}
($2 == "ns.dead_leafn") {printf("%f	%s\n",$1*0.015,$2); a=1;}
($2 == "ns.leafn_store") {printf("%f	%s\n",$1*0.015,$2); a=1;}
($2 == "ns.leafn_transfer") {printf("%f	%s\n",$1*0.015,$2); a=1;}
($2 == "ns.live_stemn") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "ns.livestemn_store") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "ns.livestemn_transfer") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "ns.dead_stemn") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "ns.deadstemn_store") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "ns.deadstemn_transfer") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "ns.live_crootn") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "ns.livecrootn_store") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "ns.livecrootn_transfer") {printf("%f	%s\n",$1*0.00,$2); a=1;}
($2 == "ns.dead_crootn") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "ns.deadcrootn_store") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "ns.deadcrootn_transfer") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "ns.frootn") {printf("%f	%s\n",$1*0.015,$2); a=1;}
($2 == "ns.frootn_store") {printf("%f	%s\n",$1*0.015,$2); a=1;}
($2 == "ns.frootn_transfer") {printf("%f	%s\n",$1*0.015,$2); a=1;}
($2 == "ns.npool") {printf("%f	%s\n",0.003,$2); a=1;}
($2 == "ns.cwdn") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "ns.retransn") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "epv.wstress_days") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "epv.max_fparabs") {printf("%f  %s\n" ,-9999,$2); a=1;}
($2 == "epv.min_vwc") {printf("%f  %s\n" ,-9999,$2); a=1;}


($2 == "veg_parm_ID") {if ($1 == "3") printf("%d  %s\n" ,2,$2);
                      else printf("%d  %s\n" ,$1,$2); a=1;}

(a == 0) {printf("%s	%s\n",$1,$2);}
