;===================================================================================================
; read in atmospheric data from CASCADE files
;
; can read in from multiple time steps and place in same array
; copy coordinate variables from CASCADE array to finished array
;
; performs temporal interpolation if flag_interp = .TRUE.
;===============================================================

function read_cascade(generic_path:string,timestep_path:integer,timestep_idx:integer,filename:string,\
                            var_no:string,variable_name:string,flag_interp:logical,t_step_new_ec:integer)

local generic_path,timestep_path,timestep_idx,filename,var_no,variable_name,flag_interp,\
      t_step_new_ec,no_time_slices,time_in,t_count,data_path_matr,data_path,in_file,\
      in_data_tmp,time_data_tmp,no_time_in,time_in,data_size,in_data_ts,var_dims,\
      data_size_new,in_data,i_dim,no_steps_tot,no_steps_per_old_step,time_out,\
      i_time,time_tmp,i_st,i_en,time_dim,in_data_ti


begin
;  setfileoption("nc","Format","NetCDF4Classic")
;  setfileoption("nc","Compression",1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;   READ IN DATA
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; 1. work out how many time slices we have in total
  no_time_slices = dimsizes(timestep_idx)
  time_in = new(no_time_slices,integer)

  do t_count = 0,dimsizes(timestep_path)-1

     data_path_matr = (/generic_path,"/p",tostring(timestep_path(t_count)),"/v0/",var_no,"/",    \
                        filename,"p",tostring(timestep_path(t_count)),"_",var_no,".nc"/)
     data_path      = str_concat(data_path_matr)
  
  ;  Input variable dimensions: (time, latitude, longitude) ;
     in_file        = addfile(data_path,"r")
     in_data_tmp    = in_file->$variable_name$

  ;  Read in time from file
     time_data_tmp  = in_file->time
     no_time_in = dimsizes(time_data_tmp)
     ; convert CASCADE time from days to seconds
     time_in(t_count) = toint( 60*60*round(time_data_tmp(timestep_idx(t_count))*24,0) )
     data_size = dimsizes(in_data_tmp)

     if (no_time_in.gt.1) then
        ; more than one timestep in file
        ; so time must be defined co-ord variable
        ; just take desired timestep
        if (dimsizes(data_size).eq.3) then     ; time,      lat, lon
           in_data_ts = in_data_tmp(timestep_idx(t_count),:,:)
           in_data_ts!0 = "lat"
           in_data_ts!1 = "lon"
        else
           if (dimsizes(data_size).eq.4) then  ; time, lev, lat, lon
              in_data_ts = in_data_tmp(timestep_idx(t_count),:,:,:)
              in_data_ts!0 = "eta_hh"
              in_data_ts!1 = "lat"
              in_data_ts!2 = "lon"
           end if
        end if
     else
        ; time is not defined co-ord variable
        if (timestep_idx(t_count).eq.0) then
           in_data_ts = in_data_tmp
           if (dimsizes(data_size).eq.2) then     ; lat, lon
              in_data_ts!0 = "lat"
              in_data_ts!1 = "lon"
           else
              if (dimsizes(data_size).eq.3) then  ; lev, lat, lon
                 in_data_ts!0 = "eta_hh"
                 in_data_ts!1 = "lat"
                 in_data_ts!2 = "lon"
              end if
           end if
        else
           print("##############################################################")
           print("#  Error: timestep_idx should be 0 if only one time in file  #")
           print("#         check the timestamps of the data                   #")
           print("##############################################################")
        end if
     end if
     copy_VarAtts(in_data_tmp,in_data_ts)
     delete([/in_data_tmp,data_size/])

     ; set up the new empty array into which we will put our multiple time data
     if (t_count.eq.0) then
        copy_VarAtts(time_data_tmp,time_in)

        var_dims = getvardims(in_data_ts)
        data_size = dimsizes(in_data_ts)
        data_size_new = array_append_record(no_time_slices,data_size,0)
        delete([/data_size/])
        data_size = data_size_new
        delete([/data_size_new/])
        in_data = new(data_size,float)
        in_data!0    = "time"
        do i_dim = 0,dimsizes(var_dims)-1
           in_data!(i_dim+1)           = in_data_ts!(i_dim)
           in_data&$in_data!(i_dim+1)$ = in_data_ts&$in_data_ts!(i_dim)$
        end do
        delete([/var_dims,i_dim,data_size/])

        ; copy over some attributes 
;        in_data@stash_code    = in_data_ts@stash_code
;        in_data@long_name     = in_data_ts@long_name
        copy_VarAtts(in_data_ts,in_data)
     end if

     ; put time slice into in_data
     data_size = dimsizes(in_data)

     if (dimsizes(data_size).eq.3) then     ; time,      lat, lon
        in_data(t_count,:,:)   = (/in_data_ts(lat|:,lon|:)/)
     else                                   ; time, lev, lat, lon
        in_data(t_count,:,:,:) = (/in_data_ts(eta_hh|:,lat|:,lon|:)/)
     end if
     
     delete([/time_data_tmp,data_path_matr,data_path,in_file,data_size,in_data_ts/])
  
  end do

  ; hmc new
  in_data&time = time_in
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;   TEMPORAL INTERPOLATION IF REQUIRED
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  if (flag_interp) then
     ; perform linear interpolation in time
     print("    perform temporal interpolation")

     no_steps_tot = round(1+ (time_in(dimsizes(time_in)-1)-time_in(0))  /t_step_new_ec,3)
     no_steps_per_old_step = round(1+(time_in(1)-time_in(0))/t_step_new_ec,3)
     time_out = new(no_steps_tot,typeof(time_in))

     ; ensure time_in times are exactly in time_out
     do i_time=0,dimsizes(time_in)-2   ; minus 1 as start at 0, minus 1 as only interpolate within times
        time_tmp = toint( fspan(time_in(i_time),time_in(i_time+1),no_steps_per_old_step) )
        i_st =    i_time *(no_steps_per_old_step-1)
        i_en = (1+i_time)*(no_steps_per_old_step-1)
        time_out(i_st:i_en) = time_tmp
        delete([/time_tmp,i_st,i_en/])
     end do

     time_dim = 0
     in_data_ti = linint1_n_Wrap (time_in,in_data,False,time_out,0,time_dim)
     in_data_ti!0 = "time"

     ; time_out in seconds
     in_data_ti&time = time_out

     delete([/time_in,no_steps_tot,time_out,in_data/])
  else
     in_data_ti = in_data
     delete([/in_data/])
  end if

  in_data_ti&time@calendar = "gregorian"
  in_data_ti&time@standard_name = "time" ;
  in_data_ti&time@units = "seconds since 2009-04-06 00:00:00" ;

  ; delete time from attributes (not dimensions)
  delete_VarAtts(in_data_ti,"time")
  in_data_ti@center = "U.K. Met Office"

  return(in_data_ti)

end




