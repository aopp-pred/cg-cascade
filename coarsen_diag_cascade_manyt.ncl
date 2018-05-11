; coarsen_diag_cascade_manyt.ncl


load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "read_cascade.ncl"
load "ecmwf_functions.ncl"
load "cascade_pressure_ec.ncl"
load "cascade_coords.ncl"
load "reduce_area_data.ncl"
load "add_to_file.ncl"
load "init_time.ncl"

begin

  ;===================================;
  ; check input information from CLI  ;
  ;===================================;

    if (.not. isdefined("pindex1") .or. \
        .not. isdefined("pindex2") .or. \
        .not. isdefined("tindex1") .or. \
        .not. isdefined("tindex2") .or. \
        .not. isdefined("lcloud")) then
      print("usage: ncl coarsen_diag_cascade_manyt.ncl pindex1=I pindex2=J tindex1=i tindex2=j lcloud=True/False")
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

    generic_path     = "/badc/cascade/data/WarmPool-4km/xfhfc"
    timestep_path    = (/ pindex1, pindex2 /)       ; i.e.. digit from folder id  "p1" etc
    timestep_idx     = (/ tindex1, tindex2 /)       ; idx within respective path
    filename         = "CASCADE_WarmPool-4km_xfhfc_"
    no_tsteps = dimsizes(timestep_path)

    ; ecmwf data
    ecmwf_path     = "/group_workspaces/jasmin2/aopp/cg-cascade/cg07/static/ecmwf_forcing"

;;----
;;  in each section define variable number and set up full path:
;;  var_no = "409"
;;  data_path_matr = (/generic_path,"/p",tostring(timestep_path),"/",filename,"p",tostring(timestep_path),"_",var_no,".nc"/)
;;----

  ;===================================;
  ; define paths to output files      ;
  ;===================================;

    path_out     = "/group_workspaces/jasmin2/aopp/cg-cascade/cg07/data/diag_cas/raw/"          ; Output directory
    file_out_tmp = (/filename,"p",\
                     tostring(timestep_path(0)),".",tostring(timestep_idx(0)),"-p",\
                     tostring(timestep_path(no_tsteps-1)),".",tostring(timestep_idx(no_tsteps-1)),\
                     "_SCM_T639.nc"/)                                                               ; Output file name
    file_out = str_concat(file_out_tmp)
    delete([/file_out_tmp,no_tsteps/])

  ;===================================;
  ; define temporal interpolation     ;
  ;===================================;

    t_step_new_ec = 15*60         ; in seconds  : 15 minutes for T639. Ensure integer. 
    flag_interp   = True

  ;===================================;
  ; define optional smoothing         ;
  ;===================================;

  smooth_flag = False
  ;smooth_flag = False

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
;; remove forcing regions at boundaries: move inwards 0.3 deg
;    set_lat_min = -19.7
;    set_lat_max =  19.7
;    set_lon_min = 42.3
;    set_lon_max = 176.7


    flag_subset = True

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
    ;delete([/ps_tmp/]) ;save for z_sfc calc

;   unweighted average because *input* fields equally spaced,
    ps_data_out = area_hi2lores_Wrap (ps_data&lon, ps_data&lat, ps_data, False,  1, lon_ec_out&lon,lat_ec_out&lat,False)

    if smooth_flag then   ;smooth data to remove unresolved features
       ps_data_out_sm = smth9_Wrap(ps_data_out,0.5,0.25,False)
       delete([/ps_data_out/])
       ps_data_out = ps_data_out_sm
       delete([/ps_data_out_sm/])
    end if

    delete([/var_no,ps_data/])


  ;=================================================================================
  ; define co-ordinate arrays

    levels =  cascade_pressure_ec(ps_data_out,ref_P0_out,\
              hyam_ec_out,hybm_ec_out,hyai_ec_out,hybi_ec_out)

    ec_pres_hybrid_levs  = levels[0]
    ec_pres_hybrid_hlevs = levels[1]
    delete([/levels/])

    time           = ec_pres_hybrid_levs&time
    time!0         = "time"
    time@long_name = "Time"
    time@calendar  = "gregorian"
    time@standard_name = "time" ;
    time@units     = "seconds since 2009-04-06 00:00:00" ;

    lat        = ec_pres_hybrid_levs&lat
    lat@units  = "deg N"
    lon        = ec_pres_hybrid_levs&lon
    lon@units  = "deg E"

  ;=================================================================;
  ; SET UP SAVING DATA TO NCL FILE                                  ;

    print("=========================================")
    print(" **  OPEN NCL FILE AND DEFINE CO-ORDS ** ")
    print("=========================================")

  ;===================================================================
  ; Define dimensions of variables
  ;

    n_lat   = dimsizes(lat)
    n_lon   = dimsizes(lon)

    system("/bin/rm -f " + path_out + file_out)    ; remove if exists
    fout  = addfile (path_out + file_out, "c")  ; open output file - create

;    fout  = addfile (path_out + file_out, "w")  ; open output file - read and write

  ;===================================================================
  ; explicitly declare file definition mode. Improve efficiency.

    setfileoption(fout,"DefineMode",True)

  ;===================================================================
  ; create global attributes of the file

    fAtt               = True            ; assign file attributes
    fAtt@title         = "CASCADE diagnostics on IFS grid: T639"
    fAtt@Conventions   = "None"
    fAtt@creation_date = systemfunc ("date")
    fileattdef( fout, fAtt )            ; copy file attributes

  ;===================================================================
  ; predefine the coordinate variables and their dimensionality
  ; Note: to get an UNLIMITED record dimension, we set the dimensionality
  ; to -1 (or the actual size) and set the dimension name to True.

    dimNames = (/"time", "lat", "lon"/)    ;;   , "nlev", "nlevp1", "nlevs"/)
    dimSizes = (/ -1   ,  n_lat,  n_lon/)  ;;   , n_lev, n_levp1, n_levs/)
    dimUnlim = (/ True ,  False,  False/)  ;;   , False,   False,  False/)
    filedimdef(fout,dimNames,dimSizes,dimUnlim)

  ;===================================================================
  ; start writing to file

    add_to_file(fout,time                 ,"time")
    add_to_file(fout,lat                  ,"lat")
    add_to_file(fout,lon                  ,"lon")

    ;=============================================================
    delete([/file_out,path_out,n_lat,n_lon/])   ;;,n_levp1,n_levs,nlev/])
    delete([/dimNames,dimSizes,dimUnlim/])
    delete([/fAtt/])

  ;==================================================================================
  ; param 1. total precipitation rate
  ;         n.b. we will take the 60-min averaged fields
  ;         though the 15-min accomulated are also available.
  ;          Cascade data in kg/m2/s

    print("total precipitation rate")

    var_no = "5216"
    variable_name = "precipitation_flux"

    ppt_tmp = read_cascade(generic_path,timestep_path,timestep_idx,filename,var_no,variable_name,flag_interp,time_out)
    ppt_data = reduce_area_data(ppt_tmp,flag_subset,set_lat_min,set_lat_max,set_lon_min,set_lon_max,3)
    delete([/ppt_tmp/])

;   unweighted average because *input* fields equally spaced,
;   NOTE CONSERVATIVE FOR PRECIPITATION
;    opt = True
;    opt@NLATo = 640
;    ppt_data_out = area_conserve_remap_Wrap (ppt_data&lon, ppt_data&lat, ppt_data, lon_ec_out&lon,lat_ec_out&lat,opt)
    ppt_data_out = area_hi2lores_Wrap (ppt_data&lon, ppt_data&lat, ppt_data, False,  1, lon_ec_out&lon,lat_ec_out&lat,False)



    add_to_file(fout,ppt_data_out         ,"precipitation_flux")

    if smooth_flag then   ;smooth data to remove unresolved features
       ppt_data_out_sm = smth9_Wrap(ppt_data_out,0.5,0.25,False)
       delete([/ppt_data_out/])
       ppt_data_out = ppt_data_out_sm
       delete([/ppt_data_out_sm/])
       add_to_file(fout,ppt_data_out         ,"precipitation_flux_sm")
    end if

    delete([/var_no,ppt_data,ppt_data_out/])

;  ;==================================================================================
;  ; param 2. outgoing longwave radiation at TOA
;  ;          Cascade data in W/m2
;
;    print("toa OLR")
;
;    var_no = "2205"
;    variable_name = "toa_outgoing_longwave_flux"
;
;    olr_tmp = read_cascade(generic_path,timestep_path,timestep_idx,filename,var_no,variable_name,flag_interp,time_out)
;    olr_data = reduce_area_data(olr_tmp,flag_subset,set_lat_min,set_lat_max,set_lon_min,set_lon_max,3)
;    delete([/olr_tmp/])
;
;;   unweighted average because *input* fields equally spaced,
;    olr_data_out = area_hi2lores_Wrap (olr_data&lon, olr_data&lat, olr_data, False,  1, lon_ec_out&lon,lat_ec_out&lat,False)
;
;    add_to_file(fout,olr_data_out         ,"toa_outgoing_longwave_flux")
;
;    if smooth_flag then   ;smooth data to remove unresolved features
;       olr_data_out_sm = smth9_Wrap(olr_data_out,0.5,0.25,False)
;       delete([/olr_data_out/])
;       olr_data_out = olr_data_out_sm
;       delete([/olr_data_out_sm/])
;       add_to_file(fout,olr_data_out         ,"toa_outgoing_longwave_flux_sm")
;    end if
;
;    delete([/var_no,olr_data,olr_data_out/])
;
  ;==================================================================================
  ; add new parameters here
  ;



    print("diagnostics extraction successfully completed")

  end


