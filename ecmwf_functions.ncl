;===================================================================================================
; read in boundary condition data from ECMWF files
;
; reads in entire file and then restricts to chosen times
; interpolates to required frequency from 6 or 12 hourly data
;         -> 1 hourly if flag_interp = .FALSE.
;         -> finer    if flag_interp = .TRUE.
;===============================================================

;============
;  CONTAINS
;  temp_interp
;  read_ecmwf
;  read_ecmwf_pad
;  blend_data


function temp_interp(in_data:float,t_required:double)

local    in_data,t_required,time_dim

begin

 ;======== OLD ===================
 ;   time_in = in_data&time
 ;   no_steps_tot = round(1+ (time_in(dimsizes(time_in)-1)-time_in(0))  /dt_new,3)
 ;   no_steps_per_old_step = round(1+(time_in(1)-time_in(0))/dt_new,3)
 ;   
 ;   time_out = new(no_steps_tot,typeof(time_in))
 ;   ; ensure time_in times are exactly in time_out
 ;   do i_time=0,dimsizes(time_in)-2   ; minus 1 as start at 0, minus 1 as only interpolate within times
 ;      time_tmp = fspan(time_in(i_time),time_in(i_time+1),no_steps_per_old_step)
 ;      i_st =    i_time *(no_steps_per_old_step-1)
 ;      i_en = (1+i_time)*(no_steps_per_old_step-1)
 ;      time_out(i_st:i_en) = time_tmp
 ;      delete([/time_tmp,i_st,i_en/])
 ;   end do
 ;
 ;   time_dim = 0
 ;   in_data_ti = linint1_n_Wrap (in_data&time,in_data,False,time_out,0,time_dim)
 ;   in_data_ti!0 = "time"
 ;
 ;   delete([/time_in,no_steps_tot,no_steps_per_old_step,time_out/])
 ;=================================

    time_dim = 0
    in_data_ti = linint1_n_Wrap (in_data&time,in_data,False,t_required,0,time_dim)
    in_data_ti!0 = "time"

    delete([/time_dim/])

    return(in_data_ti)

end


function read_ecmwf(ecmwf_path:string,lat_ec_out:float,lon_ec_out:float,\
                    flag_subset:logical,set_lat_min:float,set_lat_max:float,\
                    set_lon_min:float,set_lon_max:float,\
                    flag_interp:logical,t_required:double,t_step_new_ec:float)

local    ecmwf_path,lat_ec_out,lon_ec_out,flag_subset,set_lat_min,set_lat_max,\
         set_lon_min,set_lon_max,flag_interp,t_required,t_step_new_ec,\
         lat_min,lat_max,lon_min,lon_max,data_path_matr,data_path,\
         ec_file,lat_ec_in,lon_ec_in,time,data_tmp,data_tmp2,\
         dt_new,time_zero,dsize,t_soil_tmp,q_soil_tmp,\
         albedo,anor,heat_rough,high_veg_cover,high_veg_type, isor,\
         low_veg_cover,low_veg_type,lsm,mom_rough,open_sst,sdfor,sdor,\
         sea_ice_frct,slor,t_sea_ice,z_sfc,orog,t_skin,q_skin,t_soil,q_soil
 

begin

    ;================================================;
    ;  Read in ECMWF file for reference co-ordinates ;
    ;================================================;
    print("read in ECMWF forcing files")

    if (flag_subset) then 
       lat_min = set_lat_min 
       lat_max = set_lat_max 
       lon_min = set_lon_min 
       lon_max = set_lon_max 
    else 
       lat_min = lat_ec_out(0)
       lat_max = lat_ec_out(dimsizes(lat_ec_out)-1) 
       lon_min = lon_ec_out(0)
       lon_max = lon_ec_out(dimsizes(lon_ec_out)-1)
    end if

    time_zero = (158*365+51*366+31+28+31+5)*24

    if (flag_interp) then
       dt_new =  t_step_new_ec
    else
       dt_new = 1.0
    end if

    ;================================================;
    ;  albedo

    print("    albedo")

    data_path_matr = (/ecmwf_path,"/albedo_0001.grib"/)
    data_path      = str_concat(data_path_matr)

    ec_file = addfile(data_path,"r")

    albedo       = ec_file->AL_GDS4_SFC             ; albedo
    lat_ec_in    = ec_file->g4_lat_1                ; lat
    lon_ec_in    = ec_file->g4_lon_2                ; lon
    time         = ec_file->initial_time0_hours     ; hours since 1800-01-01 00:00

    lat_ec_in!0 = "lat"
    lon_ec_in!0 = "lon"

    ; change time so in **seconds** after UTC 6 April 2009 0000
    time = (time - time_zero)*60*60
 
    albedo!0         = "time"
    albedo&time      = time
    albedo!1         = "lat"
    albedo&lat       = lat_ec_in
    albedo!2         = "lon"
    albedo&lon       = lon_ec_in

    ; reduce area data and temporal interpolation
    data_tmp = albedo(time|:,{lat|lat_min:lat_max},{lon|lon_min:lon_max})
    delete([/albedo/])

    albedo = temp_interp(data_tmp,t_required)

    delete([/data_tmp,ec_file,lat_ec_in,lon_ec_in,time/])

    ;================================================;
    ;  anor

    print("    anor - orography AN")

    data_path_matr = (/ecmwf_path,"/anor_0001.grib"/)
    data_path      = str_concat(data_path_matr)

    ec_file = addfile(data_path,"r")

    anor         = ec_file->ANOR_GDS4_SFC             ; anor
    lat_ec_in    = ec_file->g4_lat_1                ; lat
    lon_ec_in    = ec_file->g4_lon_2                ; lon
    time         = ec_file->initial_time0_hours     ; hours since 1800-01-01 00:00

    lat_ec_in!0 = "lat"
    lon_ec_in!0 = "lon"

    ; change time so in **seconds** after UTC 6 April 2009 0000
    time = (time - time_zero)*60*60

    anor!0         = "time"
    anor&time      = time
    anor!1         = "lat"
    anor&lat  = lat_ec_in
    anor!2         = "lon"
    anor&lon = lon_ec_in

    ; reduce area data and temporal interpolation
    data_tmp = anor(time|:,{lat|lat_min:lat_max},{lon|lon_min:lon_max})
    delete([/anor/])

    anor = temp_interp(data_tmp,t_required)

    delete([/data_tmp,ec_file,lat_ec_in,lon_ec_in,time/])

    ;================================================;
    ;  heat_rough

    print("    heat_rough")

    data_path_matr = (/ecmwf_path,"/heat_rough_log_fc_0001.grib"/)
    data_path      = str_concat(data_path_matr)

    ec_file = addfile(data_path,"r")

    heat_rough   = ec_file->FLSR_GDS4_SFC             ; heat_rough
    lat_ec_in    = ec_file->g4_lat_1                ; lat
    lon_ec_in    = ec_file->g4_lon_2                ; lon
    time         = ec_file->initial_time0_hours     ; hours since 1800-01-01 00:00

    lat_ec_in!0 = "lat"
    lon_ec_in!0 = "lon"

    ; change time so in **seconds** after UTC 6 April 2009 0000
    time = (time - time_zero)*60*60

    heat_rough!0         = "time"
    heat_rough&time      = time
    heat_rough!1         = "lat"
    heat_rough&lat  = lat_ec_in
    heat_rough!2         = "lon"
    heat_rough&lon = lon_ec_in

    ; reduce area data and temporal interpolation
    data_tmp = heat_rough(time|:,{lat|lat_min:lat_max},{lon|lon_min:lon_max})
    delete([/heat_rough/])
    
    ; reading in log of heat rough - raise to exponential
    data_tmp2 = exp(data_tmp)
    copy_VarCoords(data_tmp,data_tmp2)

    heat_rough = temp_interp(data_tmp2,t_required)
    copy_VarAtts(data_tmp,heat_rough)

    delete([/data_tmp,ec_file,lat_ec_in,lon_ec_in,time/])


    ;================================================;
    ;  high_veg_cover

    print("    high_veg_cover")

    data_path_matr = (/ecmwf_path,"/high_veg_cover_0001.grib"/)
    data_path      = str_concat(data_path_matr)

    ec_file = addfile(data_path,"r")

    high_veg_cover = ec_file->CVH_GDS4_SFC             ; high_veg_cover
    lat_ec_in      = ec_file->g4_lat_1                ; lat
    lon_ec_in      = ec_file->g4_lon_2                ; lon
    time           = ec_file->initial_time0_hours     ; hours since 1800-01-01 00:00

    lat_ec_in!0 = "lat"
    lon_ec_in!0 = "lon"

    ; change time so in **seconds** after UTC 6 April 2009 0000
    time = (time - time_zero)*60*60

    high_veg_cover!0         = "time"
    high_veg_cover&time      = time
    high_veg_cover!1         = "lat"
    high_veg_cover&lat  = lat_ec_in
    high_veg_cover!2         = "lon"
    high_veg_cover&lon = lon_ec_in

    ; reduce area data and temporal interpolation
    data_tmp = high_veg_cover(time|:,{lat|lat_min:lat_max},{lon|lon_min:lon_max})
    delete([/high_veg_cover/])

    high_veg_cover = temp_interp(data_tmp,t_required)

    delete([/data_tmp,ec_file,lat_ec_in,lon_ec_in,time/])

    ;================================================;
    ;  high_veg_type

    print("    high_veg_type")

    data_path_matr = (/ecmwf_path,"/high_veg_type_0001.grib"/)
    data_path      = str_concat(data_path_matr)

    ec_file = addfile(data_path,"r")

    high_veg_type  = ec_file->TVH_GDS4_SFC             ; high_veg_type
    lat_ec_in      = ec_file->g4_lat_1                ; lat
    lon_ec_in      = ec_file->g4_lon_2                ; lon
    time           = ec_file->initial_time0_hours     ; hours since 1800-01-01 00:00

    lat_ec_in!0 = "lat"
    lon_ec_in!0 = "lon"

    ; change time so in **seconds** after UTC 6 April 2009 0000
    time = (time - time_zero)*60*60

    high_veg_type!0         = "time"
    high_veg_type&time      = time
    high_veg_type!1         = "lat"
    high_veg_type&lat  = lat_ec_in
    high_veg_type!2         = "lon"
    high_veg_type&lon = lon_ec_in

    ; reduce area data and temporal interpolation
    data_tmp = high_veg_type(time|:,{lat|lat_min:lat_max},{lon|lon_min:lon_max})
    delete([/high_veg_type/])

    high_veg_type = temp_interp(data_tmp,t_required)

    delete([/data_tmp,ec_file,lat_ec_in,lon_ec_in,time/])

    ;================================================;
    ;  isor

    print("    isor - orography IS")

    data_path_matr = (/ecmwf_path,"/isor_0001.grib"/)
    data_path      = str_concat(data_path_matr)

    ec_file = addfile(data_path,"r")

    isor           = ec_file->ISOR_GDS4_SFC           ; isor
    lat_ec_in      = ec_file->g4_lat_1                ; lat
    lon_ec_in      = ec_file->g4_lon_2                ; lon
    time           = ec_file->initial_time0_hours     ; hours since 1800-01-01 00:00

    lat_ec_in!0 = "lat"
    lon_ec_in!0 = "lon"

    ; change time so in **seconds** after UTC 6 April 2009 0000
    time = (time - time_zero)*60*60

    isor!0         = "time"
    isor&time      = time
    isor!1         = "lat"
    isor&lat  = lat_ec_in
    isor!2         = "lon"
    isor&lon = lon_ec_in

    ; reduce area data and temporal interpolation
    data_tmp = isor(time|:,{lat|lat_min:lat_max},{lon|lon_min:lon_max})
    delete([/isor/])

    isor = temp_interp(data_tmp,t_required)

    delete([/data_tmp,ec_file,lat_ec_in,lon_ec_in,time/])


    ;================================================;
    ;  low_veg_cover

    print("    low_veg_cover")

    data_path_matr = (/ecmwf_path,"/low_veg_cover_0001.grib"/)
    data_path      = str_concat(data_path_matr)

    ec_file = addfile(data_path,"r")

    low_veg_cover  = ec_file->CVL_GDS4_SFC            ; low_veg_cover
    lat_ec_in      = ec_file->g4_lat_1                ; lat
    lon_ec_in      = ec_file->g4_lon_2                ; lon
    time           = ec_file->initial_time0_hours     ; hours since 1800-01-01 00:00

    lat_ec_in!0 = "lat"
    lon_ec_in!0 = "lon"

    ; change time so in **seconds** after UTC 6 April 2009 0000
    time = (time - time_zero)*60*60

    low_veg_cover!0         = "time" 
    low_veg_cover&time      = time
    low_veg_cover!1         = "lat"
    low_veg_cover&lat  = lat_ec_in
    low_veg_cover!2         = "lon"
    low_veg_cover&lon = lon_ec_in

    ; reduce area data and temporal interpolation
    data_tmp = low_veg_cover(time|:,{lat|lat_min:lat_max},{lon|lon_min:lon_max})
    delete([/low_veg_cover/])

    low_veg_cover = temp_interp(data_tmp,t_required)

    delete([/data_tmp,ec_file,lat_ec_in,lon_ec_in,time/])

    ;================================================;
    ;  low_veg_type

    print("    low_veg_type")

    data_path_matr = (/ecmwf_path,"/low_veg_type_0001.grib"/)
    data_path      = str_concat(data_path_matr)

    ec_file = addfile(data_path,"r")

    low_veg_type   = ec_file->TVL_GDS4_SFC             ; low_veg_type
    lat_ec_in      = ec_file->g4_lat_1                ; lat
    lon_ec_in      = ec_file->g4_lon_2                ; lon
    time           = ec_file->initial_time0_hours     ; hours since 1800-01-01 00:00

    lat_ec_in!0 = "lat"
    lon_ec_in!0 = "lon"

    ; change time so in **seconds** after UTC 6 April 2009 0000
    time = (time - time_zero)*60*60

    low_veg_type!0         = "time"
    low_veg_type&time      = time
    low_veg_type!1         = "lat"
    low_veg_type&lat  = lat_ec_in
    low_veg_type!2         = "lon"
    low_veg_type&lon = lon_ec_in

    ; reduce area data and temporal interpolation
    data_tmp = low_veg_type(time|:,{lat|lat_min:lat_max},{lon|lon_min:lon_max})
    delete([/low_veg_type/])

    low_veg_type = temp_interp(data_tmp,t_required)

    delete([/data_tmp,ec_file,lat_ec_in,lon_ec_in,time/])

    ;================================================;
    ;  lsm

    print("    lsm")

    data_path_matr = (/ecmwf_path,"/lsm_0001.grib"/)
    data_path      = str_concat(data_path_matr)

    ec_file = addfile(data_path,"r")
    
    lsm            = ec_file->LSM_GDS4_SFC             ; lsm
    lat_ec_in      = ec_file->g4_lat_1                ; lat
    lon_ec_in      = ec_file->g4_lon_2                ; lon
    time           = ec_file->initial_time0_hours     ; hours since 1800-01-01 00:00

    lat_ec_in!0 = "lat"
    lon_ec_in!0 = "lon"

    ; change time so in **seconds** after UTC 6 April 2009 0000
    time = (time - time_zero)*60*60

    lsm!0         = "time"
    lsm&time      = time
    lsm!1         = "lat"
    lsm&lat  = lat_ec_in
    lsm!2         = "lon"
    lsm&lon = lon_ec_in

    ; reduce area data and temporal interpolation
    data_tmp = lsm(time|:,{lat|lat_min:lat_max},{lon|lon_min:lon_max})
    delete([/lsm/])

    lsm = temp_interp(data_tmp,t_required)

    delete([/data_tmp,ec_file,lat_ec_in,lon_ec_in,time/])

    ;================================================;
    ;  mom_rough

    print("    mom_rough")

    data_path_matr = (/ecmwf_path,"/mom_rough_fc_0001.grib"/)
    data_path      = str_concat(data_path_matr)

    ec_file = addfile(data_path,"r")
    
    mom_rough      = ec_file->FSR_GDS4_SFC             ; mom_rough
    lat_ec_in      = ec_file->g4_lat_1                ; lat
    lon_ec_in      = ec_file->g4_lon_2                ; lon
    time           = ec_file->initial_time0_hours     ; hours since 1800-01-01 00:00

    lat_ec_in!0 = "lat"
    lon_ec_in!0 = "lon"

    ; change time so in **seconds** after UTC 6 April 2009 0000
    time = (time - time_zero)*60*60

    mom_rough!0         = "time"
    mom_rough&time      = time
    mom_rough!1         = "lat"
    mom_rough&lat  = lat_ec_in
    mom_rough!2         = "lon"
    mom_rough&lon = lon_ec_in

    ; reduce area data and temporal interpolation
    data_tmp = mom_rough(time|:,{lat|lat_min:lat_max},{lon|lon_min:lon_max})
    delete([/mom_rough/])

    mom_rough = temp_interp(data_tmp,t_required)

    delete([/data_tmp,ec_file,lat_ec_in,lon_ec_in,time/])

    ;================================================;
    ;  open_sst

    print("    open_sst")

    data_path_matr = (/ecmwf_path,"/open_sst_0001.grib"/)
    data_path      = str_concat(data_path_matr)

    ec_file = addfile(data_path,"r")
    
    open_sst       = ec_file->SSTK_GDS4_SFC           ; open_sst
    lat_ec_in      = ec_file->g4_lat_1                ; lat
    lon_ec_in      = ec_file->g4_lon_2                ; lon
    time           = ec_file->initial_time0_hours     ; hours since 1800-01-01 00:00

    lat_ec_in!0 = "lat"
    lon_ec_in!0 = "lon"

    ; change time so in **seconds** after UTC 6 April 2009 0000
    time = (time - time_zero)*60*60

    open_sst!0         = "time"
    open_sst&time      = time
    open_sst!1         = "lat"
    open_sst&lat  = lat_ec_in
    open_sst!2         = "lon"
    open_sst&lon = lon_ec_in

    ; reduce area data and temporal interpolation
    data_tmp = open_sst(time|:,{lat|lat_min:lat_max},{lon|lon_min:lon_max})
    delete([/open_sst/])

    open_sst = temp_interp(data_tmp,t_required)

    delete([/data_tmp,ec_file,lat_ec_in,lon_ec_in,time/])

    ;================================================;
    ;  sdfor

    print("    sdfor - orography filtered SD")

    data_path_matr = (/ecmwf_path,"/sdfor_0001.grib"/)
    data_path      = str_concat(data_path_matr)

    ec_file = addfile(data_path,"r")
    
    sdfor          = ec_file->SDFOR_GDS4_SFC           ; sdfor
    lat_ec_in      = ec_file->g4_lat_1                ; lat
    lon_ec_in      = ec_file->g4_lon_2                ; lon
    time           = ec_file->initial_time0_hours     ; hours since 1800-01-01 00:00

    lat_ec_in!0 = "lat"
    lon_ec_in!0 = "lon"

    ; change time so in **seconds** after UTC 6 April 2009 0000
    time = (time - time_zero)*60*60

    sdfor!0         = "time"
    sdfor&time      = time
    sdfor!1         = "lat"
    sdfor&lat  = lat_ec_in
    sdfor!2         = "lon"
    sdfor&lon = lon_ec_in

    ; reduce area data and temporal interpolation
    data_tmp = sdfor(time|:,{lat|lat_min:lat_max},{lon|lon_min:lon_max})
    delete([/sdfor/])

    sdfor = temp_interp(data_tmp,t_required)

    delete([/data_tmp,ec_file,lat_ec_in,lon_ec_in,time/])

    ;================================================;
    ;  sdor

    print("    sdor - orography SD")

    data_path_matr = (/ecmwf_path,"/sdor_0001.grib"/)
    data_path      = str_concat(data_path_matr)

    ec_file = addfile(data_path,"r")
    
    sdor           = ec_file->SDOR_GDS4_SFC           ; sdor
    lat_ec_in      = ec_file->g4_lat_1                ; lat
    lon_ec_in      = ec_file->g4_lon_2                ; lon
    time           = ec_file->initial_time0_hours     ; hours since 1800-01-01 00:00

    lat_ec_in!0 = "lat"
    lon_ec_in!0 = "lon"

    ; change time so in **seconds** after UTC 6 April 2009 0000
    time = (time - time_zero)*60*60

    sdor!0         = "time"
    sdor&time      = time
    sdor!1         = "lat"
    sdor&lat  = lat_ec_in
    sdor!2         = "lon"
    sdor&lon = lon_ec_in

    ; reduce area data and temporal interpolation
    data_tmp = sdor(time|:,{lat|lat_min:lat_max},{lon|lon_min:lon_max})
    delete([/sdor/])

    sdor = temp_interp(data_tmp,t_required)

    delete([/data_tmp,ec_file,lat_ec_in,lon_ec_in,time/])

    ;================================================;
    ;  sea_ice_frct

    print("    sea_ice_frct")

    data_path_matr = (/ecmwf_path,"/sea_ice_frct_0001.grib"/)
    data_path      = str_concat(data_path_matr)

    ec_file = addfile(data_path,"r")
    
    sea_ice_frct   = ec_file->CI_GDS4_SFC             ; sea_ice_frct
    lat_ec_in      = ec_file->g4_lat_1                ; lat
    lon_ec_in      = ec_file->g4_lon_2                ; lon
    time           = ec_file->initial_time0_hours     ; hours since 1800-01-01 00:00

    lat_ec_in!0 = "lat"
    lon_ec_in!0 = "lon"

    ; change time so in **seconds** after UTC 6 April 2009 0000
    time = (time - time_zero)*60*60

    sea_ice_frct!0         = "time"
    sea_ice_frct&time      = time
    sea_ice_frct!1         = "lat"
    sea_ice_frct&lat  = lat_ec_in
    sea_ice_frct!2         = "lon"
    sea_ice_frct&lon = lon_ec_in

    ; reduce area data and temporal interpolation
    data_tmp = sea_ice_frct(time|:,{lat|lat_min:lat_max},{lon|lon_min:lon_max})
    delete([/sea_ice_frct/])

    sea_ice_frct = temp_interp(data_tmp,t_required)

    delete([/data_tmp,ec_file,lat_ec_in,lon_ec_in,time/])


    ;================================================;
    ;  slor

    print("    slor - orography SL")

    data_path_matr = (/ecmwf_path,"/slor_0001.grib"/)
    data_path      = str_concat(data_path_matr)

    ec_file = addfile(data_path,"r")
    
    slor           = ec_file->SLOR_GDS4_SFC           ; slor
    lat_ec_in      = ec_file->g4_lat_1                ; lat
    lon_ec_in      = ec_file->g4_lon_2                ; lon
    time           = ec_file->initial_time0_hours     ; hours since 1800-01-01 00:00

    lat_ec_in!0 = "lat"
    lon_ec_in!0 = "lon"

    ; change time so in **seconds** after UTC 6 April 2009 0000
    time = (time - time_zero)*60*60

    slor!0         = "time"
    slor&time      = time
    slor!1         = "lat"
    slor&lat  = lat_ec_in
    slor!2         = "lon"
    slor&lon = lon_ec_in

    ; reduce area data and temporal interpolation
    data_tmp = slor(time|:,{lat|lat_min:lat_max},{lon|lon_min:lon_max})
    delete([/slor/])

    slor = temp_interp(data_tmp,t_required)

    delete([/data_tmp,ec_file,lat_ec_in,lon_ec_in,time/])

    ;================================================;
    ;  z_sfc

    print("    z_sfc / orog")

    data_path_matr = (/ecmwf_path,"/z_sfc_0001.grib"/)
    data_path      = str_concat(data_path_matr)

    ec_file = addfile(data_path,"r")

    z_sfc          = ec_file->Z_GDS4_SFC              ; z_sfc
    lat_ec_in      = ec_file->g4_lat_1                ; lat
    lon_ec_in      = ec_file->g4_lon_2                ; lon
    time           = ec_file->initial_time0_hours     ; hours since 1800-01-01 00:00

    lat_ec_in!0 = "lat"
    lon_ec_in!0 = "lon"

    ; change time so in **seconds** after UTC 6 April 2009 0000
    time = (time - time_zero)*60*60

    z_sfc!0         = "time"
    z_sfc&time      = time
    z_sfc!1         = "lat"
    z_sfc&lat  = lat_ec_in
    z_sfc!2         = "lon"
    z_sfc&lon = lon_ec_in

    ; reduce area data and temporal interpolation
    data_tmp = z_sfc(time|:,{lat|lat_min:lat_max},{lon|lon_min:lon_max})
    delete([/z_sfc/])

    z_sfc = temp_interp(data_tmp,t_required)

    ; calculate orography from surface geopotential
    orog = z_sfc/9.80665
    copy_VarCoords(z_sfc,orog)
    copy_VarAtts(z_sfc,orog)


    delete([/data_tmp,ec_file,lat_ec_in,lon_ec_in,time/])

    ;================================================;
    ;  t_skin
 
    print("    t_skin")

    data_path_matr = (/ecmwf_path,"/t_skin_0001.grib"/)
    data_path      = str_concat(data_path_matr)

    ec_file = addfile(data_path,"r")

    t_skin       = ec_file->SKT_GDS4_SFC            ; t_skin
    lat_ec_in    = ec_file->g4_lat_1                ; lat
    lon_ec_in    = ec_file->g4_lon_2                ; lon
    time         = ec_file->initial_time0_hours     ; hours since 1800-01-01 00:00

    lat_ec_in!0 = "lat"
    lon_ec_in!0 = "lon"

    ; change time so in **seconds** after UTC 6 April 2009 0000
    time = (time - time_zero)*60*60

    t_skin!0         = "time"
    t_skin&time      = time
    t_skin!1         = "lat"
    t_skin&lat       = lat_ec_in
    t_skin!2         = "lon"
    t_skin&lon       = lon_ec_in

    ; reduce area data and temporal interpolation
    data_tmp = t_skin(time|:,{lat|lat_min:lat_max},{lon|lon_min:lon_max})
    delete([/t_skin/])

    t_skin = temp_interp(data_tmp,t_required)

    delete([/data_tmp,ec_file,lat_ec_in,lon_ec_in,time/])

    ;================================================;
    ;  q_skin

    print("    q_skin")
 
    data_path_matr = (/ecmwf_path,"/q_skin_0001.grib"/)
    data_path      = str_concat(data_path_matr)
 
    ec_file = addfile(data_path,"r")

    q_skin       = ec_file->SRC_GDS4_SFC             ; q_skin
    lat_ec_in    = ec_file->g4_lat_1                ; lat
    lon_ec_in    = ec_file->g4_lon_2                ; lon 
    time         = ec_file->initial_time0_hours     ; hours since 1800-01-01 00:00
 
    lat_ec_in!0 = "lat" 
    lon_ec_in!0 = "lon" 
 
    ; change time so in **seconds** after UTC 6 April 2009 0000
    time = (time - time_zero)*60*60 
 
    q_skin!0         = "time"
    q_skin&time      = time
    q_skin!1         = "lat"
    q_skin&lat       = lat_ec_in
    q_skin!2         = "lon"
    q_skin&lon       = lon_ec_in

    ; reduce area data and temporal interpolation
    q_skin_full = q_skin
    data_tmp = q_skin(time|:,{lat|lat_min:lat_max},{lon|lon_min:lon_max})
    delete([/q_skin/])

    q_skin = temp_interp(data_tmp,t_required)

    delete([/data_tmp,ec_file,lat_ec_in,lon_ec_in,time/])

    ;================================================;
    ;  t_sea_ice

    print("    t_sea_ice")

    dsize = dimsizes(q_skin_full)
    t_sea_ice_tmp   = new((/dsize(0),dsize(1),dsize(2),4/),typeof(t_skin))

    data_path_matr = (/ecmwf_path,"/t_sea_ice_0001.grib"/)
    data_path      = str_concat(data_path_matr)
    ec_file        = addfile(data_path,"r")
    t_sea_ice_tmp(:,:,:,0)    = ec_file->ISTL1_GDS4_DBLY         ; t_sea_ice - layer 1
    t_sea_ice_tmp(:,:,:,1)    = ec_file->ISTL2_GDS4_DBLY         ; t_sea_ice - layer 2
    t_sea_ice_tmp(:,:,:,2)    = ec_file->ISTL3_GDS4_DBLY         ; t_sea_ice - layer 3
    t_sea_ice_tmp(:,:,:,3)    = ec_file->ISTL4_GDS4_DBLY         ; t_sea_ice - layer 4
    lat_ec_in                 = ec_file->g4_lat_1                ; lat
    lon_ec_in                 = ec_file->g4_lon_2                ; lon
    time                      = ec_file->initial_time0_hours     ; hours since 1800-01-01 00:00

    lat_ec_in!0 = "lat"
    lon_ec_in!0 = "lon"

    ; change time so in **seconds** after UTC 6 April 2009 0000
    time = (time - time_zero)*60*60

    t_sea_ice_tmp!0         = "time"
    t_sea_ice_tmp&time      = time
    t_sea_ice_tmp!1         = "lat"
    t_sea_ice_tmp&lat       = lat_ec_in
    t_sea_ice_tmp!2         = "lon"
    t_sea_ice_tmp&lon       = lon_ec_in
    t_sea_ice_tmp!3         = "nlevs"
    t_sea_ice_tmp&nlevs     = ispan(1,4,1)

    t_sea_ice = t_sea_ice_tmp(time | :, nlevs | :, lat | :, lon | :)

    ; reduce area data and temporal interpolation
    data_tmp = t_sea_ice(time|:,nlevs|:,{lat|lat_min:lat_max},{lon|lon_min:lon_max})
    delete([/t_sea_ice/])

    t_sea_ice = temp_interp(data_tmp,t_required)

    delete([/data_tmp,ec_file,lat_ec_in,lon_ec_in,time/])
    delete([/t_sea_ice_tmp/])


    ;================================================;
    ;  t_soil
 
    print("    t_soil")

    dsize = dimsizes(q_skin_full)
    delete([/q_skin_full/])
    t_soil_tmp   = new((/dsize(0),dsize(1),dsize(2),4/),typeof(t_skin))
 
    data_path_matr = (/ecmwf_path,"/t_soil_1_0001.grib"/)
    data_path      = str_concat(data_path_matr)
    ec_file        = addfile(data_path,"r")
    t_soil_tmp(:,:,:,0)   = ec_file->STL1_GDS4_DBLY          ; t_soil
    lat_ec_in             = ec_file->g4_lat_1                ; lat
    lon_ec_in             = ec_file->g4_lon_2                ; lon
    time                  = ec_file->initial_time0_hours     ; hours since 1800-01-01 00:00
    delete([/data_path_matr,data_path,ec_file/])

    data_path_matr = (/ecmwf_path,"/t_soil_2_0001.grib"/)
    data_path      = str_concat(data_path_matr)
    ec_file        = addfile(data_path,"r")
    t_soil_tmp(:,:,:,1)   = ec_file->STL2_GDS4_DBLY          ; t_soil
    delete([/data_path_matr,data_path,ec_file/])

    data_path_matr = (/ecmwf_path,"/t_soil_3_0001.grib"/)
    data_path      = str_concat(data_path_matr)
    ec_file        = addfile(data_path,"r")
    t_soil_tmp(:,:,:,2)   = ec_file->STL3_GDS4_DBLY          ; t_soil
    delete([/data_path_matr,data_path,ec_file/])

    data_path_matr = (/ecmwf_path,"/t_soil_4_0001.grib"/)
    data_path      = str_concat(data_path_matr)
    ec_file        = addfile(data_path,"r")
    t_soil_tmp(:,:,:,3)   = ec_file->STL4_GDS4_DBLY         ; t_soil
    delete([/data_path_matr,data_path,ec_file/])

    lat_ec_in!0 = "lat"
    lon_ec_in!0 = "lon"
 
    ; change time so in **seconds** after UTC 6 April 2009 0000
    time = (time - time_zero)*60*60
 
    t_soil_tmp!0         = "time"
    t_soil_tmp&time      = time
    t_soil_tmp!1         = "lat"
    t_soil_tmp&lat       = lat_ec_in
    t_soil_tmp!2         = "lon"
    t_soil_tmp&lon       = lon_ec_in
    t_soil_tmp!3         = "nlevs"
    t_soil_tmp&nlevs     = ispan(1,4,1)
 
    t_soil = t_soil_tmp(time | :, nlevs | :, lat | :, lon | :)

    ; reduce area data and temporal interpolation
    data_tmp = t_soil(time|:,nlevs|:,{lat|lat_min:lat_max},{lon|lon_min:lon_max})
    delete([/t_soil/])
 
    t_soil = temp_interp(data_tmp,t_required)
 
    delete([/data_tmp,lat_ec_in,lon_ec_in,time/])
    delete([/t_soil_tmp/])

    ;================================================;
    ;  q_soil

    print("    q_soil")

    q_soil_tmp   = new((/dsize(0),dsize(1),dsize(2),4/),typeof(q_skin))

    data_path_matr = (/ecmwf_path,"/q_soil_1_0001.grib"/)
    data_path      = str_concat(data_path_matr)
    ec_file        = addfile(data_path,"r")
    q_soil_tmp(:,:,:,0)   = ec_file->SWVL1_GDS4_DBLY         ; q_soil
    lat_ec_in             = ec_file->g4_lat_1                ; lat
    lon_ec_in             = ec_file->g4_lon_2                ; lon
    time                  = ec_file->initial_time0_hours     ; hours since 1800-01-01 00:00
    delete([/data_path_matr,data_path,ec_file/])

    data_path_matr = (/ecmwf_path,"/q_soil_2_0001.grib"/)
    data_path      = str_concat(data_path_matr)
    ec_file        = addfile(data_path,"r")
    q_soil_tmp(:,:,:,1)   = ec_file->SWVL2_GDS4_DBLY         ; q_soil
    delete([/data_path_matr,data_path,ec_file/])

    data_path_matr = (/ecmwf_path,"/q_soil_3_0001.grib"/)
    data_path      = str_concat(data_path_matr)
    ec_file        = addfile(data_path,"r")
    q_soil_tmp(:,:,:,2)   = ec_file->SWVL3_GDS4_DBLY         ; q_soil
    delete([/data_path_matr,data_path,ec_file/])

    data_path_matr = (/ecmwf_path,"/q_soil_4_0001.grib"/)
    data_path      = str_concat(data_path_matr)
    ec_file        = addfile(data_path,"r")
    q_soil_tmp(:,:,:,3)   = ec_file->SWVL4_GDS4_DBLY         ; q_soil
    delete([/data_path_matr,data_path,ec_file/])

    lat_ec_in!0 = "lat"
    lon_ec_in!0 = "lon"

    ; change time so in **seconds** after UTC 6 April 2009 0000
    time = (time - time_zero)*60*60

    q_soil_tmp!0         = "time"
    q_soil_tmp&time      = time
    q_soil_tmp!1         = "lat"
    q_soil_tmp&lat       = lat_ec_in
    q_soil_tmp!2         = "lon"
    q_soil_tmp&lon       = lon_ec_in
    q_soil_tmp!3         = "nlevs"
    q_soil_tmp&nlevs     = ispan(1,4,1)

    q_soil = q_soil_tmp(time | :, nlevs | :, lat | :, lon | :)

    ; reduce area data and temporal interpolation
    data_tmp = q_soil(time|:,nlevs|:,{lat|lat_min:lat_max},{lon|lon_min:lon_max})
    delete([/q_soil/])

    q_soil = temp_interp(data_tmp,t_required)

    delete([/data_tmp,lat_ec_in,lon_ec_in,time/])
    delete([/q_soil_tmp,dsize/])


    return([/ albedo, anor, heat_rough, high_veg_cover, high_veg_type, isor, low_veg_cover, low_veg_type, lsm, mom_rough, open_sst, sdfor, sdor, sea_ice_frct, slor, t_sea_ice, z_sfc, orog, t_skin, q_skin, t_soil, q_soil /])

end

;===================================================================================================
; read in ERAI 3D fields to pad out upper layers of CASCADE data
;
;===============================================================

function read_ecmwf_pad(ecmwf_path:string,ecmwf_var:string,lat_ec_out:float,lon_ec_out:float,\
                    flag_subset:logical,set_lat_min:float,set_lat_max:float,\
                    set_lon_min:float,set_lon_max:float,\
                    flag_interp:logical,t_required:double,t_step_new_ec:float)

local    ecmwf_path,ecmwf_var,lat_ec_out,lon_ec_out,flag_subset,set_lat_min,set_lat_max,\
         set_lon_min,set_lon_max,flag_interp,t_required,t_step_new_ec,\
         lat_min,lat_max,lon_min,lon_max,data_path_matr,data_path,\
         ec_file,lat_ec_in,lon_ec_in,time,hyam,hybm,data_tmp,\
         time_zero,erai_data


begin

    ;================================================;
    ;  Read in ECMWF file for reference co-ordinates ;
    ;================================================;
    print("     read in ERAi files for padding")

    if (flag_subset) then
       lat_min = set_lat_min
       lat_max = set_lat_max
       lon_min = set_lon_min
       lon_max = set_lon_max
    else
       lat_min = lat_ec_out(0)
       lat_max = lat_ec_out(dimsizes(lat_ec_out)-1)
       lon_min = lon_ec_out(0)
       lon_max = lon_ec_out(dimsizes(lon_ec_out)-1)
    end if

    time_zero = (158*365+51*366+31+28+31+5)*24

    ; read in file
    data_path_matr = (/ecmwf_path,"/",ecmwf_var,"_ml_0001_20090406-20090421_LA.grib"/)
    data_path      = str_concat(data_path_matr)

    ec_file = addfile(data_path,"r")

    if (ecmwf_var.eq."U") then
       erai_data    = ec_file->U_GDS4_HYBL          ; erai_data = U
    else
       if (ecmwf_var.eq."V") then
          erai_data    = ec_file->V_GDS4_HYBL          ; erai_data = V
       else
          if (ecmwf_var.eq."T") then
             erai_data    = ec_file->T_GDS4_HYBL          ; erai_data = T
          else
             if (ecmwf_var.eq."W") then
                erai_data    = ec_file->W_GDS4_HYBL          ; erai_data = W
             else
              if (ecmwf_var.eq."Q") then
                 erai_data    = ec_file->Q_GDS4_HYBL          ; erai_data = Q
                 erai_data = 0.*erai_data                     ; want array of zeros
              else
                if (ecmwf_var.eq."QL") then
                   erai_data    = ec_file->CLWC_GDS4_HYBL      ; erai_data = QL                ; want array of zeros
                   erai_data = 0.*erai_data                    ; want array of zeros
                else
                   if (ecmwf_var.eq."QI") then
                      erai_data    = ec_file->CIWC_GDS4_HYBL      ; erai_data = QI                ; want array of zeros
                      erai_data = 0.*erai_data                    ; want array of zeros
                   else
                      if (ecmwf_var.eq."CF") then
                         erai_data    = ec_file->CC_GDS4_HYBL      ; erai_data = CF
                         erai_data = 0.*erai_data                  ; want array of zeros
                      else
                      end if
                   end if
                end if
              end if
             end if
          end if
       end if
    end if

    lat_ec_in    = ec_file->g4_lat_4                ; lat
    lon_ec_in    = ec_file->g4_lon_5                ; lon
    time         = ec_file->initial_time0_hours     ; hours since 1800-01-01 00:00
    hyam         = ec_file->lv_HYBL1_a
    hybm         = ec_file->lv_HYBL1_b

    lat_ec_in!0 = "lat"
    lon_ec_in!0 = "lon"
    hyam!0      = "nlev"
    hybm!0      = "nlev"

    ; change time so in **seconds** after UTC 6 April 2009 0000
    time = (time - time_zero)*60*60

    erai_data!0         = "time"
    erai_data&time      = time
    erai_data!1         = "nlev"
    erai_data&nlev      = hyam&nlev
    erai_data!2         = "lat"
    erai_data&lat       = lat_ec_in
    erai_data!3         = "lon"
    erai_data&lon       = lon_ec_in

    ; reduce area data and temporal interpolation
    data_tmp = erai_data(time|:,nlev|:,{lat|lat_min:lat_max},{lon|lon_min:lon_max})
    delete([/erai_data/])

    erai_data = temp_interp(data_tmp,t_required)

    delete([/data_tmp,ec_file,lat_ec_in,lon_ec_in,time/])

    return(erai_data)

end


;===================================================================================================
; blends data from CASCADE and ERAI to create coherent profile
;===============================================================

function blend_data(in_data:double,in_data_erai:float,no_blend_over:float)

local    in_data,in_data_erai,no_blend_over,dsize,frac_data,frac_erai,missingvals,d_blend,i_blend,\
         frac_data_4d,frac_erai_4d,blended_data,in1D,indices

begin

     dsize = dimsizes(in_data)

     frac_data = new(dsize(1),typeof(in_data))
     frac_erai = new(dsize(1),typeof(in_data))

     missingvals = ismissing(in_data)

     d_blend = 1/no_blend_over 
     i_blend = 1

     do i_lev = 0,dsize(1)-1
       if any(missingvals(:,i_lev,:,:)) then ; data missing in CASCADE in_data
         frac_data(i_lev) = 0.0
         frac_erai(i_lev) = 1.0
       else
         if (i_blend.le.no_blend_over) then
            frac_data(i_lev) = i_blend*d_blend
            frac_erai(i_lev) = 1.0 - i_blend*d_blend
            i_blend = i_blend+1
         else
            frac_data(i_lev) = 1.0
            frac_erai(i_lev) = 0.0
         end if
       end if
     end do

     ; set fill values to zero
     in1D = ndtooned(in_data)
     indices = ind(ismissing(in1D))
     in1D(indices) = 0
     in_data = onedtond(in1D,dimsizes(in_data))
     delete([/in1D,indices/])

     ; blend using fractions
     frac_data_4d = conform_dims(dsize,frac_data,1)
     frac_erai_4d = conform_dims(dsize,frac_erai,1)

     blended_data = in_data*frac_data_4d + in_data_erai*frac_erai_4d

    delete([/dsize,frac_data,frac_erai,missingvals,d_blend,i_blend/])
    delete([/frac_data_4d,frac_erai_4d/])

    return(blended_data)

end




