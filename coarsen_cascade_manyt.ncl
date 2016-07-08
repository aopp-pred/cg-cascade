; coarsen_cascade_manyt.ncl


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "read_cascade.ncl"
load "calc_advection.ncl"
load "ecmwf_functions.ncl"
load "cascade_coords.ncl"
load "reduce_area_data.ncl"
load "cascade_pressure_ec.ncl"
load "add_to_file.ncl"
load "cascade_etadotdpdeta.ncl"
load "cascade_ttend.ncl"
load "cascade_qtend.ncl"
load "cascade_utend.ncl"
load "cascade_vtend.ncl"
load "cascade_geopotential.ncl"
load "cascade_geostrophic.ncl"
load "init_time.ncl"
load "read_cascade_missingcloud.ncl"

begin

  ;===================================;
  ; check input information from CLI  ;
  ;===================================;

    if (.not. isdefined("pindex1") .or. \
        .not. isdefined("pindex2") .or. \
        .not. isdefined("tindex1") .or. \
        .not. isdefined("tindex2") .or. \
        .not. isdefined("lcloud")) then
      print("usage: ncl coarsen_cascade_manyt.ncl pindex1=I pindex2=J tindex1=i tindex2=j lcloud=True/False")
      exit()
    end if

    if ((pindex1 .lt. 1) .or. (pindex1 .gt. 120)) then
      print("error: pindex1 must be in the range: 1 <= pindex1 <= 120")
      exit()
    end if
    if ((pindex2 .lt. 1) .or. (pindex2 .gt. 120)) then
      print("error: pindex2 must be in the range: 1 <= pindex2 <= 120")
      exit()
    end if
    if ((tindex1 .ne. 0) .and. (tindex1 .ne. 1)) then
      print("error: tindex1 must be 0 or 1")
      exit()
    end if
    if ((tindex2 .ne. 0) .and. (tindex2 .ne. 1)) then
      print("error: tindex2 must be 0 or 1")
      exit()
    end if

    if (pindex1 .eq. pindex2) then
      if ((tindex1 .ne. 0) .or. (tindex2 .ne.1)) then
        print("when pindex1=pindex2 require tindex1=0, tindex2=1")
        exit()
      end if
    end if
    if (pindex1 .ne. pindex2) then
      if (pindex2 .ne. (pindex1 + 1)) then
        print("when pindex1/=pindex2 require pindex2=pindex1+1")
        exit()
      end if
      if ((tindex1 .ne. 1) .or. (tindex2 .ne. 0)) then
        print("when pindex2=pindex1+1 require tindex1=1, tindex2=0")
        exit()
      end if
    end if

    if (.not.(typeof(lcloud).eq."logical"))
      print("lcloud must be logical: True or False")
      exit()
    end if


  ;===================================;
  ; define paths to input variables   ;
  ;===================================;

;    generic_path     = "/badc/cascade/data/WarmPool-4km/xfhfc" ;JASMIN
    generic_path     = "../DATA/WarmPool-4km/FULLDATA" ;CIRRUS
    timestep_path    = (/ pindex1, pindex2 /)       ; i.e.. digit from folder id  "p1" etc
    timestep_idx     = (/ tindex1, tindex2 /)       ; idx within respective path
    filename         = "CASCADE_WarmPool-4km_xfhfc_"
    no_tsteps = dimsizes(timestep_path)

    ; ecmwf data
    ecmwf_path     = "../DATA/ForcingFilesT639"

;;----
;;  in each section define variable number and set up full path:
;;  var_no = "409"
;;  data_path_matr = (/generic_path,"/p",tostring(timestep_path),"/",filename,"p",tostring(timestep_path),"_",var_no,".nc"/)
;;----

  ;===================================;
  ; define paths to output files      ;
  ;===================================;

;    path_out     = "/group_workspaces/jasmin2/aopp/christensen/make_scm_input/DATA/SCM_IN_T639/"          ; Output directory JASMIN
    path_out     = "../DATA/SCM_IN_T639/"          ; Output directory CIRRUS
    file_out_tmp = (/filename,"p",\
                     tostring(timestep_path(0)),".",tostring(timestep_idx(0)),"-p",\
                     tostring(timestep_path(no_tsteps-1)),".",tostring(timestep_idx(no_tsteps-1)),\
                     "_test_new_time_SCM_T639.nc"/)                                                               ; Output file name
    file_out = str_concat(file_out_tmp)
    delete([/file_out_tmp,no_tsteps/])

  ;===================================;
  ; define temporal interpolation     ;
  ;===================================;

    t_step_new_ec = 15*60         ; in seconds  : 15 minutes for T639. Ensure integer. 
    flag_interp   = True

  ;===================================;
  ; select subset data for testing    ;
  ;===================================;

;; reduced?
;    set_lat_min = -5.0
;    set_lat_max =  5.0
;    set_lon_min = 50.5
;    set_lon_max = 64.5
;; full?
    set_lat_min = -20.0
    set_lat_max =  20.0
    set_lon_min = 42.0
    set_lon_max = 177.0
;-- orography?
;    set_lat_min = -20.0
;    set_lat_max = -10.0
;    set_lon_min = 45.0
;    set_lon_max = 55.0


    flag_subset = True

  ;===================================================;
  ;  determine cloud frac timesteps from timestep_idx ;
  ;===================================================;

    dsize = dimsizes(timestep_path)
    no_files = dsize(0)
    idx_useme = new(dsize,integer)
    do i=0,no_files-1
       idx_useme(i) = 1
    end do
    delete([/i/])
    flag_offset = False
    do idx = 0,no_files-1
       if (timestep_path(idx).eq.1).and.(timestep_idx(idx).eq.1)
          print("== Cloud file offset from others ==")
          idx_useme(idx) = 0
          flag_offset = True
       end if 
    end do
    if (.not.flag_offset)
       if (timestep_path(0).ge.2) ; still need to create offset
          print("== Cloud file offset from others ==")
          idx_useme(0) = 0
          flag_offset = True
       end if
    end if
    timestep_path_cf = timestep_path(ind(idx_useme.eq.1))
    timestep_idx_cf  = timestep_idx(ind(idx_useme.eq.1))
 
    if flag_offset   ; have deleted one file from input, add new on end
       timestep_path_cf_tmp = new(dsize,integer)
       timestep_idx_cf_tmp  = new(dsize,integer)
       timestep_path_cf_tmp(0:no_files-2) = timestep_path_cf
       timestep_idx_cf_tmp(0:no_files-2)  = timestep_idx_cf
 
       ; add extras on end
       last_idx     = timestep_idx(no_files-1)
       last_pathno  = timestep_path(no_files-1)
       if (last_idx.eq.0) 
          timestep_path_cf_tmp(no_files-1) = timestep_path(no_files-1)
          timestep_idx_cf_tmp(no_files-1)  = 1
       else
          timestep_path_cf_tmp(no_files-1) = timestep_path(no_files-1)+1
          timestep_idx_cf_tmp(no_files-1)  = 0
       end if
       delete([/timestep_path_cf,timestep_idx_cf/])
       timestep_path_cf = timestep_path_cf_tmp
       timestep_idx_cf  = timestep_idx_cf_tmp
       delete([/timestep_path_cf_tmp,timestep_idx_cf_tmp/])
    end if
  ;================================================;
  ;  Read in ECMWF file for reference co-ordinates ;
  ;================================================;

    ecmwf_coords = cascade_coords(flag_subset,set_lat_min,set_lat_max,set_lon_min,set_lon_max)
    hyam_ec_out    = ecmwf_coords[0]
    hybm_ec_out    = ecmwf_coords[1]
    hyai_ec_out    = ecmwf_coords[2]
    hybi_ec_out    = ecmwf_coords[3]
    lat_ec_out     = ecmwf_coords[4]
    lon_ec_out     = ecmwf_coords[5]
    ref_P0_out     = ecmwf_coords[6]
    delete([/ecmwf_coords/])

  ;==========================
  ;  SET UP INTERPOLATION
  ;==========================

  ;==================================================================================
  ; param -1. determine output time vector from trusted files

    print("determine out time")

    var_no = "409"
    variable_name = "surface_air_pressure"

    time_out = init_time(generic_path,timestep_path,timestep_idx,filename,var_no,variable_name,flag_interp,t_step_new_ec)


  ;==================================================================================
  ; param 0. surface pressure after timestep
  ;         (going to need this for inter-
  ;          polation to hybrid levels)
  ;          Cascade data in Pa

    print("surface pressure calculation")

    var_no = "409"
    variable_name = "surface_air_pressure"

    ps_tmp = read_cascade(generic_path,timestep_path,timestep_idx,filename,var_no,variable_name,flag_interp,time_out)
    ps_data = reduce_area_data(ps_tmp,flag_subset,set_lat_min,set_lat_max,set_lon_min,set_lon_max,3)
    delete([/ps_tmp/])

;   unweighted average because *input* fields equally spaced,
    ps_data_out = area_hi2lores_Wrap (ps_data&lon, ps_data&lat, ps_data, False,  1, lon_ec_out&lon,lat_ec_out&lat,False)

    delete([/var_no,ps_data/])

  ;=================================================================================
  ;  define pressure on ECMWF full and half levels

    levels =  cascade_pressure_ec(ps_data_out,ref_P0_out,\
              hyam_ec_out,hybm_ec_out,hyai_ec_out,hybi_ec_out)

    ec_pres_hybrid_levs  = levels[0]
    ec_pres_hybrid_hlevs = levels[1]
    delete([/levels/])

  ; define co-ordinate arrays
    time           = ec_pres_hybrid_levs&time
    time!0         = "time"
    time@long_name = "Time"
    time@calendar  = "gregorian"
    time@standard_name = "time" ;
    time@units     = "seconds since 2009-04-06 00:00:00" ;
    nlev           = ec_pres_hybrid_levs&nlev
    nlev@units     = "count"
    nlev@long_name = "Atmospheric Model Levels"
    nlevp1         = hyai_ec_out&nlevp1
    nlevp1@units   = "count"
    nlevp1@long_name = "Atmospheric Model Half Levels"
    nlevs          = ispan(1,4,1)
    nlevs!0        = "nlevs"
    nlevs@units    = "count"
    nlevs@long_name = "Soil/Sea-Ice Model Levels"

    lat        = ec_pres_hybrid_levs&lat
    lat@units  = "deg N"
    lon        = ec_pres_hybrid_levs&lon
    lon@units  = "deg E"

    ;;-- set up initial time for input file
    ; date and second => initial time
    ; work out what the date is
    ; time is in seconds since 2009-04-06 00:00:00
;
;    start_date   = 20090406
;    start_seconds = 0
;    date = time+start_date               ;floor(time/(60*60*24))+start_date
;    date@long_name = "Date"
;    date@calendar  = "gregorian"
;    date@standard_name = "date" ;
;    date@units     = "s" ;
;    second = date*0+start_seconds
;    second@long_name     = "Second"
;    second@calendar      = "gregorian"
;    second@standard_name = "second" ;
;    second@units         = "seconds" ;
;    date!0   = "time"
;    second!0 = "time"
;    delete([/start_seconds,start_date/])

  ;=================================================================;
  ; SET UP SAVING DATA TO NCL FILE                                  ;

    print("=========================================")
    print(" **  OPEN NCL FILE AND DEFINE CO-ORDS ** ")
    print("=========================================")

  ;===================================================================
  ; Define dimensions of variables
  ;

    n_lev   = dimsizes(nlev)
    n_lat   = dimsizes(lat)
    n_lon   = dimsizes(lon)
    n_levp1 = dimsizes(nlevp1)
    n_levs  = dimsizes(nlevs)

    system("/bin/rm -f " + path_out + file_out)    ; remove if exists
    fout  = addfile (path_out + file_out, "c")  ; open output file - create

;    fout  = addfile (path_out + file_out, "w")  ; open output file - read and write

  ;===================================================================
  ; explicitly declare file definition mode. Improve efficiency.

    setfileoption(fout,"DefineMode",True)

  ;===================================================================
  ; create global attributes of the file

    fAtt               = True            ; assign file attributes
    fAtt@title         = "SCM input file: T639"
    fAtt@Conventions   = "None"
    fAtt@creation_date = systemfunc ("date")
    fileattdef( fout, fAtt )            ; copy file attributes

  ;===================================================================
  ; predefine the coordinate variables and their dimensionality
  ; Note: to get an UNLIMITED record dimension, we set the dimensionality
  ; to -1 (or the actual size) and set the dimension name to True.

    dimNames = (/"time", "lat", "lon", "nlev", "nlevp1", "nlevs"/)
    dimSizes = (/ -1   ,  n_lat,  n_lon, n_lev, n_levp1, n_levs/)
    dimUnlim = (/ True ,  False,  False, False,   False,  False/)
    filedimdef(fout,dimNames,dimSizes,dimUnlim)

  ;===================================================================
  ; start writing to file

    add_to_file(fout,time                 ,"time")
    add_to_file(fout,nlev                 ,"nlev")
    add_to_file(fout,nlevp1               ,"nlevp1")
    add_to_file(fout,nlevs                ,"nlevs")
    add_to_file(fout,lat                  ,"lat")
    add_to_file(fout,lon                  ,"lon")
    add_to_file(fout,ec_pres_hybrid_levs  , "pressure_f")
    add_to_file(fout,ec_pres_hybrid_hlevs , "pressure_h")
    add_to_file(fout,ps_data_out          , "ps")
;    add_to_file(fout,date                 , "date")
;    add_to_file(fout,second               , "second")

    ;=============================================================
    delete([/file_out,path_out,n_lat,n_lon,n_levp1,n_levs,nlev/])
    delete([/dimNames,dimSizes,dimUnlim/])
;    delete([/second,date/])
    delete([/ps_data_out,fAtt,nlevs/])

  ;==================================================================================
  ; param 1. pressure on rho levels after timestep
  ;         (going to need this for inter-
  ;          polation to hybrid levels)  
  ;          note: ordering is surface to top of atmos

    print("rho level pressure calculation")

    var_no = "407"
    variable_name = "air_pressure"
    prho_tmp = read_cascade(generic_path,timestep_path,timestep_idx,filename,var_no,variable_name,flag_interp,time_out)
    prho_data = reduce_area_data(prho_tmp,flag_subset,set_lat_min,set_lat_max,set_lon_min,set_lon_max,4)
    delete([/prho_tmp/])

    prho_data_out = area_hi2lores_Wrap (prho_data&lon, prho_data&lat, prho_data, False,  1, lon_ec_out&lon,lat_ec_out&lat, False)

    delete([/var_no,variable_name,prho_data/])

  ;==================================================================================
  ; param 2. pressure on theta levels after timestep
  ;         (going to need this for inter-
  ;          polation to hybrid levels)
  ;          note: ordering is surface to top of atmos

    print("theta level pressure calculation")
  
    var_no = "408"  
    variable_name = "air_pressure"
    ptheta_tmp = read_cascade(generic_path,timestep_path,timestep_idx,filename,var_no,variable_name,flag_interp,time_out)
    ptheta_data = reduce_area_data(ptheta_tmp,flag_subset,set_lat_min,set_lat_max,set_lon_min,set_lon_max,4)
    delete([/ptheta_tmp/])

    ptheta_data_out = area_hi2lores_Wrap (ptheta_data&lon, ptheta_data&lat, ptheta_data, False,  1, lon_ec_out&lon,lat_ec_out&lat, False)

    delete([/var_no,variable_name/])

  ;=========================================
  ;  INPUT ECMWF BOUNDARY FILES
  ;
  ;    -> feed in required timesteps time_out

    ecmwf_data = read_ecmwf(ecmwf_path,lat_ec_out,lon_ec_out,flag_subset,set_lat_min,set_lat_max,\
                      set_lon_min,set_lon_max,time_out)

    albedo          = ecmwf_data[0]
    anor            = ecmwf_data[1]
    heat_rough      = ecmwf_data[2]
    high_veg_cover  = ecmwf_data[3]
    high_veg_type   = ecmwf_data[4]
    isor            = ecmwf_data[5]
    low_veg_cover   = ecmwf_data[6]
    low_veg_type    = ecmwf_data[7]
    lsm             = ecmwf_data[8]
    mom_rough       = ecmwf_data[9]
    open_sst        = ecmwf_data[10]
    sdfor           = ecmwf_data[11]
    sdor            = ecmwf_data[12]
    sea_ice_frct    = ecmwf_data[13]
    slor            = ecmwf_data[14]
    t_sea_ice       = ecmwf_data[15]
    z_sfc           = ecmwf_data[16]
    orog            = ecmwf_data[17]
    t_skin          = ecmwf_data[18]
    q_skin          = ecmwf_data[19]
    t_soil          = ecmwf_data[20]
    q_soil          = ecmwf_data[21]

    delete([/ecmwf_data/])

    ; create snow variable
    snow = 0.0*t_skin
    copy_VarCoords(t_skin,snow)
    copy_VarAtts(t_skin,snow)
    snow@long_name = "Snow Depth"
    snow@units     = "m"

    ; if over land point, set open_sst to t_skin
    if(any(ismissing(open_sst))) then
       open_sst_1D = ndtooned(open_sst)
       t_skin_1D   = ndtooned(t_skin)
       open_sst_1D(ind(ismissing(open_sst_1D))) = t_skin_1D(ind(ismissing(open_sst_1D)))
       open_sst = onedtond(open_sst_1D,dimsizes(open_sst))
       delete([/open_sst_1D,t_skin_1D/])
    end if

    ;-- save to file
    add_to_file(fout,albedo         ,"albedo")
    add_to_file(fout,anor           ,"anor")
    add_to_file(fout,heat_rough     ,"heat_rough")
    add_to_file(fout,high_veg_cover ,"high_veg_cover")
    add_to_file(fout,high_veg_type  ,"high_veg_type")
    add_to_file(fout,isor           ,"isor")
    add_to_file(fout,low_veg_cover  ,"low_veg_cover")
    add_to_file(fout,low_veg_type   ,"low_veg_type")
    add_to_file(fout,lsm            ,"lsm")
    add_to_file(fout,mom_rough      ,"mom_rough")
    add_to_file(fout,open_sst       ,"open_sst")
    add_to_file(fout,sdfor          ,"sdfor")
    add_to_file(fout,sdor           ,"sdor")
    add_to_file(fout,sea_ice_frct   ,"sea_ice_frct")
    add_to_file(fout,slor           ,"slor")
    add_to_file(fout,t_sea_ice      ,"t_sea_ice")
    add_to_file(fout,snow           ,"snow")
    add_to_file(fout,orog           ,"orog")
;    add_to_file(fout,t_skin         ,"t_skin")  ;HMC don't write to file: now reading in from CASCADE
    add_to_file(fout,q_skin         ,"q_skin")
    add_to_file(fout,t_soil         ,"t_soil")
    add_to_file(fout,q_soil         ,"q_soil")
    add_to_file(fout,z_sfc          ,"z_sfc")

    delete([/albedo,anor,heat_rough,high_veg_cover,high_veg_type,isor,low_veg_cover/])
    delete([/low_veg_type,lsm,mom_rough,open_sst,sdfor,sdor,sea_ice_frct,slor,t_sea_ice,orog,snow/])
    delete([/t_skin,q_skin,t_soil,q_soil,z_sfc/])

  ;==================================================================================
  ; param  3. Temperature
  ;        CASCADE: theta levels
  ;        ECMWF  : full levels

    print("T calculation")

    var_no = "16004"
    variable_name = "data"
    t_tmp = read_cascade(generic_path,timestep_path,timestep_idx,filename,var_no,variable_name,flag_interp,time_out)
    t_data = reduce_area_data(t_tmp,flag_subset,set_lat_min,set_lat_max,set_lon_min,set_lon_max,4)
    delete([/t_tmp/])

    ; 1. area average
    t_data_areav = area_hi2lores_Wrap (t_data&lon, t_data&lat, t_data, False,  1, lon_ec_out&lon,lat_ec_out&lat, False)

    ; 2. vertical interpolation
    ; T is on *theta* coordinates
    linlog = 2   ; ln(p) interpolation
    lev_dim = 1
    t_data_out = int2p_n(ptheta_data_out,t_data_areav,ec_pres_hybrid_levs,linlog,lev_dim)

    t_data_out!0         = "time"
    t_data_out&time      = t_data_areav&time
    t_data_out!1         = "nlev"
    t_data_out&nlev      = ec_pres_hybrid_levs&nlev
    t_data_out!2         = "lat"
    t_data_out&lat  = t_data_areav&lat
    t_data_out!3         = "lon"
    t_data_out&lon = t_data_areav&lon

    t_data_out_erai = read_ecmwf_pad(ecmwf_path,"T",lat_ec_out,lon_ec_out,flag_subset,set_lat_min,set_lat_max,\
                      set_lon_min,set_lon_max,time_out)

    no_blend_over = 5.0
    t_data_new = blend_data(t_data_out,t_data_out_erai,no_blend_over)

    copy_VarCoords(t_data_out,t_data_new)
    copy_VarAtts(t_data_areav,t_data_new)
    delete([/t_data_out,t_data_areav/])
    ;- rename attribute 'coordinates' as clashes with Andrew Dawson scmtiles code
    t_data_new@vertical_coordinates = t_data_new@coordinates
    delete(t_data_new@coordinates)

    ;-- save to file
    add_to_file(fout,t_data_new,"t")

    delete([/var_no,variable_name,linlog,lev_dim,t_data_out_erai,no_blend_over,t_data_new/])


  ;==================================================================================
  ; param  4. Vertical velocity
  ;       - do first as need fine p field for w_to_omega
  ;        CASCADE: theta levels - in m/s upward
  ;        ECMWF  : full levels - in Pa/s

    print("W calculation")

    var_no = "150"
    variable_name = "upward_air_velocity"
    w_tmp = read_cascade(generic_path,timestep_path,timestep_idx,filename,var_no,variable_name,flag_interp,time_out)
    w_data_tmp = reduce_area_data(w_tmp,flag_subset,set_lat_min,set_lat_max,set_lon_min,set_lon_max,4)
    delete([/w_tmp/])

    ; 0. convert from m/s to Pa/s
    w_data = w_to_omega(w_data_tmp,ptheta_data,t_data)
    delete([/w_data_tmp,ptheta_data,t_data/])

    ; 1. area average
    w_data_areav = area_hi2lores_Wrap (w_data&lon, w_data&lat, w_data, False,  1, lon_ec_out&lon,lat_ec_out&lat, False)

    ; 2. vertical interpolation
    ; W is on *theta* coordinates
    linlog = 2   ; ln(p) interpolation
    lev_dim = 1
    w_data_out = int2p_n(ptheta_data_out,w_data_areav,ec_pres_hybrid_levs,linlog,lev_dim)

    w_data_out!0         = "time"
    w_data_out&time      = w_data_areav&time
    w_data_out!1         = "nlev"
    w_data_out&nlev      = ec_pres_hybrid_levs&nlev
    w_data_out!2         = "lat"
    w_data_out&lat  = w_data_areav&lat
    w_data_out!3         = "lon"
    w_data_out&lon = w_data_areav&lon

    w_data_out_erai = read_ecmwf_pad(ecmwf_path,"W",lat_ec_out,lon_ec_out,flag_subset,set_lat_min,set_lat_max,\
                      set_lon_min,set_lon_max,time_out)

    no_blend_over = 5.0
    w_data_new = blend_data(w_data_out,w_data_out_erai,no_blend_over)

    copy_VarCoords(w_data_out,w_data_new)
    copy_VarAtts(w_data_areav,w_data_new)

    delete([/var_no,variable_name,w_data,w_data_areav,linlog,lev_dim,w_data_out_erai,no_blend_over,w_data_out/])

    ;-- save to file
    add_to_file(fout,w_data_new,"omega")
    delete([/w_data_new/])

  ;==================================================================================
  ; param 5. Horizontal velocity: U
  ;        CASCADE: rho levels
  ;        ECMWF  : full levels

    print("U calculation")

    var_no = "2"
    variable_name = "eastward_wind"
    u_tmp = read_cascade(generic_path,timestep_path,timestep_idx,filename,var_no,variable_name,flag_interp,time_out)
    u_data = reduce_area_data(u_tmp,flag_subset,set_lat_min,set_lat_max,set_lon_min,set_lon_max,4)
    delete([/u_tmp/])

    ; 1. area average
    u_data_areav = area_hi2lores_Wrap (u_data&lon, u_data&lat, u_data, False,  1, lon_ec_out&lon,lat_ec_out&lat, False)

    ; 2. vertical interpolation
    ; U is on *rho* coordinates
    linlog = 2   ; ln(p) interpolation
    lev_dim = 1
    u_data_out = int2p_n(prho_data_out,u_data_areav,ec_pres_hybrid_levs,linlog,lev_dim)

    u_data_out!0         = "time"
    u_data_out&time      = u_data_areav&time
    u_data_out!1         = "nlev"
    u_data_out&nlev      = ec_pres_hybrid_levs&nlev
    u_data_out!2         = "lat"
    u_data_out&lat  = u_data_areav&lat
    u_data_out!3         = "lon"
    u_data_out&lon = u_data_areav&lon

    u_data_out_erai = read_ecmwf_pad(ecmwf_path,"U",lat_ec_out,lon_ec_out,flag_subset,set_lat_min,set_lat_max,\
                      set_lon_min,set_lon_max,time_out)

    no_blend_over = 5.0
    u_data_new = blend_data(u_data_out,u_data_out_erai,no_blend_over)
 
    copy_VarCoords(u_data_out,u_data_new)
    copy_VarAtts(u_data_areav,u_data_new)
    ;- rename attribute 'coordinates' as clashes with Andrew Dawson scmtiles code
    u_data_new@vertical_coordinates = u_data_new@coordinates
    delete(u_data_new@coordinates)

    delete([/var_no,variable_name,u_data,u_data_areav,linlog,lev_dim,u_data_out_erai,no_blend_over,u_data_out/])

   ;-- save to file
    add_to_file(fout,u_data_new,"u")
    delete([/u_data_new/])

  ;==================================================================================
  ; param 6. Horizontal velocity: V
  ;         CASCADE: rho levels
  ;         ECMWF  : full levels

    print("V calculation")

    var_no = "3"
    variable_name = "northward_wind"
    v_tmp  = read_cascade(generic_path,timestep_path,timestep_idx,filename,var_no,variable_name,flag_interp,time_out)
    v_data = reduce_area_data(v_tmp,flag_subset,set_lat_min,set_lat_max,set_lon_min,set_lon_max,4)
    delete([/v_tmp/])

    ; 1. area average
    v_data_areav = area_hi2lores_Wrap (v_data&lon, v_data&lat, v_data, False,  1, lon_ec_out&lon,lat_ec_out&lat, False)

    ; 2. vertical interpolation
    ; V is on *rho* coordinates
    linlog = 2   ; ln(p) interpolation
    lev_dim = 1
    v_data_out = int2p_n(prho_data_out,v_data_areav,ec_pres_hybrid_levs,linlog,lev_dim)

    v_data_out!0         = "time"
    v_data_out&time      = v_data_areav&time
    v_data_out!1         = "nlev"
    v_data_out&nlev      = ec_pres_hybrid_levs&nlev
    v_data_out!2         = "lat"
    v_data_out&lat  = v_data_areav&lat
    v_data_out!3         = "lon"
    v_data_out&lon = v_data_areav&lon

    v_data_out_erai = read_ecmwf_pad(ecmwf_path,"V",lat_ec_out,lon_ec_out,flag_subset,set_lat_min,set_lat_max,\
                      set_lon_min,set_lon_max,time_out)

    no_blend_over = 5.0
    v_data_new = blend_data(v_data_out,v_data_out_erai,no_blend_over)
    copy_VarCoords(v_data_out,v_data_new)
    copy_VarAtts(v_data_areav,v_data_new)
    ;- rename attribute 'coordinates' as clashes with Andrew Dawson scmtiles code
    v_data_new@vertical_coordinates = v_data_new@coordinates
    delete(v_data_new@coordinates)

    delete([/var_no,variable_name,v_data,v_data_areav,linlog,lev_dim,v_data_out_erai,no_blend_over,v_data_out/])

   ;-- save to file
    add_to_file(fout,v_data_new,"v")
    delete([/v_data_new/])

  ;==================================================================================
  ; params 7.1  Water vapour mixing ratio: q
  ;        7.2  Liquid water mixing ratio: ql
  ;        7.3  Ice water mixing ratio   : qi
  ;             CASCADE: theta levels
  ;             ECMWF  : full levels;

    print("q/ql/qi calculation")

  ;  Q: is specific humidity the same as water vapour mixing ratio? No.
  ;  CASCADE: specific humidity = mass of water vapour per unit mass of moist air (dry+vapour+liquid+ice)
  ;  IFS:     water vapour mixing ratio - mass of water vapour per unit mass of dry air

    ;== read in water vapour
    var_no = "10"
    variable_name = "specific_humidity"
    q_tmp  = read_cascade(generic_path,timestep_path,timestep_idx,filename,var_no,variable_name,flag_interp,time_out)
    q_data = reduce_area_data(q_tmp,flag_subset,set_lat_min,set_lat_max,set_lon_min,set_lon_max,4)
    delete([/q_tmp/])

    ; 1.1 area average
    q_data_areav  = area_hi2lores_Wrap (q_data&lon, q_data&lat,  q_data, False,  1, lon_ec_out&lon,lat_ec_out&lat, False)
    delete([/q_data/])

    ;== read in liquid water
    var_no = "254"
    variable_name = "mass_fraction_of_cloud_liquid_water_in_air"
    ql_tmp  = read_cascade(generic_path,timestep_path,timestep_idx,filename,var_no,variable_name,flag_interp,time_out)
    ql_data = reduce_area_data(ql_tmp,flag_subset,set_lat_min,set_lat_max,set_lon_min,set_lon_max,4)
    delete([/ql_tmp/])

    ; 1.1 area average
    ql_data_areav  = area_hi2lores_Wrap (ql_data&lon, ql_data&lat,  ql_data, False,  1, lon_ec_out&lon,lat_ec_out&lat, False)
    delete([/ql_data/])
                                                              
    ;== read in ice water
    var_no = "12"
    variable_name = "mass_fraction_of_cloud_ice_in_air"
    qi_tmp  = read_cascade(generic_path,timestep_path,timestep_idx,filename,var_no,variable_name,flag_interp,time_out)
    qi_data = reduce_area_data(qi_tmp,flag_subset,set_lat_min,set_lat_max,set_lon_min,set_lon_max,4)
    delete([/qi_tmp/])

    ; 1.1 area average
    qi_data_areav  = area_hi2lores_Wrap (qi_data&lon, qi_data&lat,  qi_data, False,  1, lon_ec_out&lon,lat_ec_out&lat, False)
    delete([/qi_data/])

    ; deduce specific mass of dry air
    qd_data_areav = 1 - (q_data_areav + ql_data_areav + qi_data_areav)

    ; calculate area av of q first, then calculate q_mr (mixingrat)
    ; to ensure conservation of moist static energy
    ; calculated difference between two routines (this vs av[calc q_mr])
    ;   - order of one thousandth of a percent in the mean,
    ;   - order 10^-4 of a percent in the max.

    ; convert to m.r.
    q_mr_data_areav  =  q_data_areav/qd_data_areav
    ql_mr_data_areav = ql_data_areav/qd_data_areav
    qi_mr_data_areav = qi_data_areav/qd_data_areav
    
    copy_VarCoords(q_data_areav, q_mr_data_areav)
    copy_VarCoords(q_data_areav,ql_mr_data_areav)
    copy_VarCoords(q_data_areav,qi_mr_data_areav)

    delete([/q_data_areav,ql_data_areav,qi_data_areav,qd_data_areav/])

    ; 1.2 vertical interpolation
    ; q is on *theta* coordinates
    linlog = 2   ; ln(p) interpolation
    lev_dim = 1
    q_mr_data_out  = int2p_n(ptheta_data_out, q_mr_data_areav,ec_pres_hybrid_levs,linlog,lev_dim)
    ql_mr_data_out = int2p_n(ptheta_data_out,ql_mr_data_areav,ec_pres_hybrid_levs,linlog,lev_dim)
    qi_mr_data_out = int2p_n(ptheta_data_out,qi_mr_data_areav,ec_pres_hybrid_levs,linlog,lev_dim)

    q_mr_data_out!0         = "time"
    q_mr_data_out&time      = q_mr_data_areav&time
    q_mr_data_out!1         = "nlev"
    q_mr_data_out&nlev      = ec_pres_hybrid_levs&nlev
    q_mr_data_out!2         = "lat"
    q_mr_data_out&lat       = q_mr_data_areav&lat
    q_mr_data_out!3         = "lon"
    q_mr_data_out&lon       = q_mr_data_areav&lon

    copy_VarCoords(q_mr_data_out,ql_mr_data_out)
    copy_VarCoords(q_mr_data_out,qi_mr_data_out)

    ; pad using erai data
    q_data_out_erai  = read_ecmwf_pad(ecmwf_path,"Q",lat_ec_out,lon_ec_out,flag_subset,set_lat_min,set_lat_max,\
                       set_lon_min,set_lon_max,time_out)
    ql_data_out_erai = read_ecmwf_pad(ecmwf_path,"QL",lat_ec_out,lon_ec_out,flag_subset,set_lat_min,set_lat_max,\
                       set_lon_min,set_lon_max,time_out)
    qi_data_out_erai = read_ecmwf_pad(ecmwf_path,"QI",lat_ec_out,lon_ec_out,flag_subset,set_lat_min,set_lat_max,\
                       set_lon_min,set_lon_max,time_out)

    no_blend_over = 5.0
    q_data_new  = blend_data(q_mr_data_out, q_data_out_erai ,no_blend_over)
    ql_data_new = blend_data(ql_mr_data_out,ql_data_out_erai,no_blend_over)
    qi_data_new = blend_data(qi_mr_data_out,qi_data_out_erai,no_blend_over)

    copy_VarCoords(q_mr_data_out, q_data_new)
    copy_VarCoords(ql_mr_data_out,ql_data_new)
    copy_VarCoords(qi_mr_data_out,qi_data_new)

    q_data_new@long_name    = "Water vapor mixing ratio"
    q_data_new@units        = "kg/kg"
    q_data_new@cell_methods = "longitude: latitude: mean"
    q_data_new@center       = "U.K. Met Office"

    copy_VarAtts(q_data_new,ql_data_new)
    ql_data_new@long_name   = "Liquid water mixing ratio"

    copy_VarAtts(q_data_new,qi_data_new)
    qi_data_new@long_name   = "Ice water mixing ratio"

    delete([/q_mr_data_out,ql_mr_data_out,qi_mr_data_out/])

    q_mr_data_out  = q_data_new
    ql_mr_data_out = ql_data_new
    qi_mr_data_out = qi_data_new

    delete([/var_no,variable_name,linlog,lev_dim/])
    delete([/q_mr_data_areav,ql_mr_data_areav,qi_mr_data_areav/])
    delete([/q_data_out_erai,ql_data_out_erai,qi_data_out_erai/])
    delete([/no_blend_over,q_data_new,ql_data_new,qi_data_new/])   

    ;-- save to file
    add_to_file(fout,q_mr_data_out   ,"q")
    add_to_file(fout,ql_mr_data_out  ,"ql")
    add_to_file(fout,qi_mr_data_out  ,"qi")

    delete([/q_mr_data_out,ql_mr_data_out,qi_mr_data_out/])

  ;==================================================================================
  ; params 8 etadotdpdeta

    cascade_etadotdpdeta(fout,hyam_ec_out,hybm_ec_out)

  ;==================================================================================
  ; param 9. Advective Tendencies
  ;         CASCADE: calculate this from the coarse grained T, U, V fields
  ;         ECMWF  : full levels

    cascade_ttend(fout)
    cascade_qtend(fout)
    cascade_utend(fout)
    cascade_vtend(fout)

  ;==================================================================================
  ; param 10. Geopotential height
 
    tv = cascade_geopotential(fout,hyam_ec_out,hybm_ec_out,hyai_ec_out,hybi_ec_out,ref_P0_out)


   ;==================================================================================
   ; param 11 "geostrophic winds" : calculate forcing term from 
   ; pressure gradient and gradient in phi
   ;

    cascade_geostrophic(fout,tv,hyam_ec_out,hybm_ec_out)

    delete([/tv/])


   ;==================================================================================
   ; param 12. Cloud fraction
   ;        CASCADE: theta levels
   ;        ECMWF  : full levels

    print("cloud fraction calculation")

   ; cloud fraction in CASCADE is split into four components:
   ; large scale ice;  large scale water
   ; convective ice ;  convective water
    ;== read in large scale liquid
    var_no = "2312"
    variable_name = "data"
    if (lcloud)
       print("lcloud True: special case for cloud files")
       cf_lsl_tmp  = read_cascade_missingcloud(generic_path,timestep_path_cf,timestep_idx_cf,filename,var_no,variable_name,flag_interp,time_out)
    else
       cf_lsl_tmp  = read_cascade(generic_path,timestep_path_cf,timestep_idx_cf,filename,var_no,variable_name,flag_interp,time_out)
    end if
    cf_lsl_data = reduce_area_data(cf_lsl_tmp,flag_subset,set_lat_min,set_lat_max,set_lon_min,set_lon_max,4)
    delete([/cf_lsl_tmp/])

    ; 1. area average
    cf_lsl_data_areav = area_hi2lores_Wrap (cf_lsl_data&lon,cf_lsl_data&lat,cf_lsl_data,False, 1,lon_ec_out&lon,lat_ec_out&lat, False)
    delete([/cf_lsl_data/])

    ;== read in large scale ice
    var_no = "2313"
    variable_name = "data"
    if (lcloud)
       print("lcloud True: special case for cloud files")
       cf_lsi_tmp  = read_cascade_missingcloud(generic_path,timestep_path_cf,timestep_idx_cf,filename,var_no,variable_name,flag_interp,time_out)
    else
       cf_lsi_tmp  = read_cascade(generic_path,timestep_path_cf,timestep_idx_cf,filename,var_no,variable_name,flag_interp,time_out)
    end if
    cf_lsi_data = reduce_area_data(cf_lsi_tmp,flag_subset,set_lat_min,set_lat_max,set_lon_min,set_lon_max,4)
    delete([/cf_lsi_tmp/])

    ; 1. area average
    cf_lsi_data_areav = area_hi2lores_Wrap (cf_lsi_data&lon,cf_lsi_data&lat,cf_lsi_data,False, 1,lon_ec_out&lon,lat_ec_out&lat, False)
    delete([/cf_lsi_data/])

    ;== read in convective liquid
    var_no = "2314"
    variable_name = "data"
    if (lcloud)
       print("lcloud True: special case for cloud files")
       cf_csl_tmp  = read_cascade_missingcloud(generic_path,timestep_path_cf,timestep_idx_cf,filename,var_no,variable_name,flag_interp,time_out)
    else
       cf_csl_tmp  = read_cascade(generic_path,timestep_path_cf,timestep_idx_cf,filename,var_no,variable_name,flag_interp,time_out)
    end if
    cf_csl_data = reduce_area_data(cf_csl_tmp,flag_subset,set_lat_min,set_lat_max,set_lon_min,set_lon_max,4)
    delete([/cf_csl_tmp/])

    ; 1. area average
    cf_csl_data_areav = area_hi2lores_Wrap (cf_csl_data&lon,cf_csl_data&lat,cf_csl_data,False,1,lon_ec_out&lon,lat_ec_out&lat, False)
    delete([/cf_csl_data/])

    ;== read in convective ice
    var_no = "2315"
    variable_name = "data"
    if (lcloud)
       print("lcloud True: special case for cloud files")
       cf_csi_tmp  = read_cascade_missingcloud(generic_path,timestep_path_cf,timestep_idx_cf,filename,var_no,variable_name,flag_interp,time_out)
    else
       cf_csi_tmp  = read_cascade(generic_path,timestep_path_cf,timestep_idx_cf,filename,var_no,variable_name,flag_interp,time_out)
    end if
    cf_csi_data = reduce_area_data(cf_csi_tmp,flag_subset,set_lat_min,set_lat_max,set_lon_min,set_lon_max,4)
    delete([/cf_csi_tmp/])

    ; 1. area average
    cf_csi_data_areav = area_hi2lores_Wrap (cf_csi_data&lon,cf_csi_data&lat,cf_csi_data,False,1,lon_ec_out&lon,lat_ec_out&lat, False)
    delete([/cf_csi_data/])

    cf_data_areav = cf_lsl_data_areav+cf_lsi_data_areav+cf_csl_data_areav+cf_csi_data_areav
    copy_VarCoords(cf_csi_data_areav,cf_data_areav)    ; time, eta_h, lat, lon

    delete([/cf_lsl_data_areav,cf_lsi_data_areav,cf_csl_data_areav,cf_csi_data_areav/])

    ; nb STRANGE OFFSET TIMESTEPS and only one timestep for first file

    ; 2. vertical interpolation
    ; T is on *theta* coordinates
    linlog = 2   ; ln(p) interpolation
    lev_dim = 1  
    ; now read in extra time steps for cf - select relevant ones
    time_list_cf = ec_pres_hybrid_levs&time
    time_st = time_list_cf(0)
    time_en = time_list_cf(dimsizes(time_list_cf)-1)

    cf_data_out = int2p_n(ptheta_data_out,cf_data_areav({time|time_st:time_en},eta_hh|:,lat|:,lon|:),ec_pres_hybrid_levs,linlog,lev_dim)

    cf_data_out!0         = "time"
    cf_data_out&time      = ec_pres_hybrid_levs&time
    cf_data_out!1         = "nlev"
    cf_data_out&nlev      = ec_pres_hybrid_levs&nlev
    cf_data_out!2         = "lat"
    cf_data_out&lat       = cf_data_areav&lat
    cf_data_out!3         = "lon"
    cf_data_out&lon       = cf_data_areav&lon

    cf_data_out_erai = read_ecmwf_pad(ecmwf_path,"CF",lat_ec_out,lon_ec_out,flag_subset,set_lat_min,set_lat_max,\
                      set_lon_min,set_lon_max,time_out)

    no_blend_over = 5.0
    cf_data_new = blend_data(cf_data_out({time|time_st:time_en},nlev|:,lat|:,lon|:),cf_data_out_erai({time|time_st:time_en},nlev|:,lat|:,lon|:),no_blend_over)

    copy_VarCoords(cf_data_out,cf_data_new)
    cf_data_new@long_name = "cloud fraction"
    cf_data_new@units     = "0-1"
    delete([/cf_data_out/])

    cf_data_out = cf_data_new

    ; set max and min to 1 and 0
    ;cf_data_1D = ndtooned(cf_data_new)
    ;print(ind((cf_data_1D.ge.1.0)))
    ;cf_data_1D(ind((cf_data_1D.ge.1.0))) = 1.0
    ;cf_data_1D(ind((cf_data_1D.le.0.0))) = 0.0
    ;cf_data_out = onedtond(cf_data_1D,dimsizes(cf_data_new))
    ;delete([/cf_data_1D/])

    delete([/var_no,variable_name,linlog,lev_dim,time_list_cf,time_st,time_en/])
    delete([/cf_data_out_erai,cf_data_areav,no_blend_over,cf_data_new/])

    add_to_file(fout, cf_data_out, "cloud_fraction")
    delete([/cf_data_out/])

  ;==================================================================================
  ; param 13. surface sensible heat flux

    print("surface sensible heat flux")

    var_no = "3217"
    variable_name = "surface_upward_sensible_heat_flux"

    sshf_tmp  = read_cascade(generic_path,timestep_path,timestep_idx,filename,var_no,variable_name,flag_interp,time_out)
    sshf_data = reduce_area_data(sshf_tmp,flag_subset,set_lat_min,set_lat_max,set_lon_min,set_lon_max,3)
    delete([/sshf_tmp/])

;   unweighted average because *input* fields equally spaced,
    sshf_data_out = area_hi2lores_Wrap (sshf_data&lon, sshf_data&lat, sshf_data, False,  1, lon_ec_out&lon,lat_ec_out&lat, False)

    ;--FLIP FLUXES--
    sshf_data_out=-1.0*sshf_data_out
    sshf_data_out@standard_name = "surface_downward_sensible_heat_flux"

    ;-- save to file
    add_to_file(fout, sshf_data_out, "sfc_sens_flx")
    delete([/sshf_data_out,var_no,sshf_data/])

  ;==================================================================================
  ; param  14. surface latent heat flux

    print("surface latent heat flux")

    var_no = "3234"
    variable_name = "surface_upward_latent_heat_flux"

    slhf_tmp  = read_cascade(generic_path,timestep_path,timestep_idx,filename,var_no,variable_name,flag_interp,time_out)
    slhf_data = reduce_area_data(slhf_tmp,flag_subset,set_lat_min,set_lat_max,set_lon_min,set_lon_max,3)
    delete([/slhf_tmp/])

;   unweighted average because *input* fields equally spaced,
    slhf_data_out = area_hi2lores_Wrap (slhf_data&lon, slhf_data&lat, slhf_data, False,  1, lon_ec_out&lon,lat_ec_out&lat, False)
    ;--FLIP FLUXES--
    slhf_data_out=-1.0*slhf_data_out
    slhf_data_out@standard_name = "surface_downward_latent_heat_flux"

    ;-- save to file
    add_to_file(fout, slhf_data_out, "sfc_lat_flx")
    delete([/slhf_data_out,var_no,slhf_data/])

  ;==================================================================================
  ; param  15. surface temperature (skin temperature)

    print("surface temperature")

    var_no = "24"
    variable_name = "surface_temperature"

    st_tmp  = read_cascade(generic_path,timestep_path,timestep_idx,filename,var_no,variable_name,flag_interp,time_out)
    st_data = reduce_area_data(st_tmp,flag_subset,set_lat_min,set_lat_max,set_lon_min,set_lon_max,3)
    delete([/st_tmp/])

;   unweighted average because *input* fields equally spaced,
    st_data_out = area_hi2lores_Wrap (st_data&lon, st_data&lat, st_data, False,  1, lon_ec_out&lon,lat_ec_out&lat, False)

    ;-- save to file
    add_to_file(fout, st_data_out, "t_skin")
    delete([/st_data_out,var_no,st_data/])



  end


