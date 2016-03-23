; cascade_coords.ncl

function cascade_coords(flag_subset:logical,set_lat_min:float,set_lat_max:float,set_lon_min:float,set_lon_max:float)

local    flag_subset, set_lat_min, set_lat_max, set_lon_min, set_lon_max, hyam_ec_out, hybm_ec_out, hyai_ec_out, hybi_ec_out, lat_ec_out, lon_ec_out, ref_P0_out

begin

  ;================================================;
  ;  Read in ECMWF file for reference co-ordinates ;
  ;================================================;
    print("read in ECMWF co-ordinates")

    ec_file = addfile("EC_ref/T_ml_0001_20090406_0000_LA.grib","r")

    hyam_ec_out  = ec_file->lv_HYBL0_a              ; \  hybrid height co-ordinate (m)
    hybm_ec_out  = ec_file->lv_HYBL0_b              ; /    at *full* pressure levels = midpoints
    hyai_ec_out  = ec_file->lv_HYBL_i1_a            ; \  hybrid height co-ordinate (m)
    hybi_ec_out  = ec_file->lv_HYBL_i1_b            ; /    at *half* pressure levels = interfaces
    lat_ec_out   = ec_file->g4_lat_3                ; lat
    lon_ec_out   = ec_file->g4_lon_4                ; lon
    ref_P0_out   = ec_file->P0                      ; reference pressure

    lat_ec_out!0 = "lat"
    lon_ec_out!0 = "lon"
    hyam_ec_out!0 = "nlev"
    hybm_ec_out!0 = "nlev"
    hyai_ec_out!0 = "nlevp1"
    hybi_ec_out!0 = "nlevp1"

    hyam_ec_out&nlev   = ispan(1,dimsizes(hyam_ec_out),1)
    hybm_ec_out&nlev   = ispan(1,dimsizes(hybm_ec_out),1)
    hyai_ec_out&nlevp1 = ispan(1,dimsizes(hyai_ec_out),1)
    hybi_ec_out&nlevp1 = ispan(1,dimsizes(hybi_ec_out),1)

    if (flag_subset) then
       print("    reduce area data")
       data_tmp = lat_ec_out({lat|set_lat_min:set_lat_max})
       delete([/lat_ec_out/])
       lat_ec_out = data_tmp
       delete([/data_tmp/])

       data_tmp = lon_ec_out({lon|set_lon_min:set_lon_max})
       delete([/lon_ec_out/])
       lon_ec_out = data_tmp
       delete([/data_tmp/])
    end if

    delete([/ec_file/])

    return([/ hyam_ec_out, hybm_ec_out, hyai_ec_out, hybi_ec_out, lat_ec_out, lon_ec_out, ref_P0_out /])

end
