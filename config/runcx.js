{
    "defaults" :
    {
        "jobname"       : "ClimberX",
        "email"         : "USER@pik-potsdam.de",
        "group"         : "megarun", 
        "omp"           : 16,
        "wall"          : 24, 
        "qos"           : "priority",
        "partition"     : "haswell",
        "job_template"  : "config/pik_submit_slurm"
    },  

    "exe_aliases" : 
        {   "climber" : "climber.x",
            "mapping" : "mapping.x"
        },
    
    "grp_aliases" :
        {   "ctl"  : "control",
            "atm"  : "atm_par",
            "bgc"  : "bgc_par",
            "geo"  : "geo_par",
            "imo"  : "imo_par",
            "lnd"  : "lnd_par",
            "ocn"  : "ocn_par",
            "sic"  : "sic_par",
            "smb"  : "smb_par",
            "hose" : "hyster_hosing"
        },

    "par_paths" :
        {   "ctl"  : "nml/control.nml",
            "atm"  : "nml/atm_par.nml",
            "lnd"  : "nml/lnd_par.nml",
            "sic"  : "nml/sic_par.nml",
            "bgc"  : "nml/bgc_par.nml",
            "geo"  : "nml/geo_par.nml",
            "imo"  : "nml/imo_par.nml",
            "ocn"  : "nml/ocn_par.nml",
            "smb"  : "nml/smb_par.nml",
            "icey" : "nml/ice_yelmo_par.nml",
            "ices" : "nml/ice_sico_par.nml",
            "hyst" : "nml/hyster_ctrl.nml"
        },

    "files" : 
        ["nml/ice_grids.nml"], 

    "dir-special" :
    {
        "restart/vilma" : "vilma_restart"
    },

    "links" : 
        ["input", "maps", "restart", "ice_data"],

    "job_queues" :
        {   "priority" :
            {   "wall" : 24  },
            "short" :
            {   "wall" : 24  },
            "medium" :
            {   "wall" : 168 },
            "long" :
            {   "wall" : 720 }
        }
}
